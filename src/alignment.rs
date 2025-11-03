use crate::io::parameters::AlignmentParameters;
use crate::models::score_matrix::MatrixType::{Ix, Iy, M};
use crate::models::score_matrix::{MatrixType, Pointer};
use crate::models::AlignGrid;
use crate::utils::{Epsilon, CHUNK_SIZE};
use num_traits::Zero;
use std::collections::HashSet;
use std::error::Error;
use std::fmt::Display;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::str::FromStr;


fn find_traceback_start<T: Copy + Display + Epsilon + FromStr + PartialOrd + Zero>(
    align_grid: &AlignGrid<T>,
    alignment_parameters: &AlignmentParameters<T>,
) -> (T, HashSet<Pointer>) {
    let mut max_val;
    let mut max_loc = HashSet::new();

    if alignment_parameters.global_alignment {
        let max_row = align_grid.m_matrix.nrow - 1;
        let max_col = align_grid.m_matrix.ncol - 1;

        let m = align_grid.m_matrix.get_score(max_row, max_col);
        max_val = m;
        max_loc.insert((M, max_row, max_col));

        let ix = align_grid.ix_matrix.get_score(max_row, max_col);
        if ix > max_val && !T::fuzzy_equals(ix, max_val) {
            max_val = ix;
            max_loc.clear();
            max_loc.insert((Ix, max_row, max_col));
        } else if T::fuzzy_equals(ix, max_val) {
            max_loc.insert((Ix, max_row, max_col));
        }

        let iy = align_grid.iy_matrix.get_score(max_row, max_col);
        if iy > max_val && !T::fuzzy_equals(iy, max_val) {
            max_val = iy;
            max_loc.clear();
            max_loc.insert((Iy, max_row, max_col));
        } else if T::fuzzy_equals(iy, max_val) {
            max_loc.insert((Iy, max_row, max_col));
        }
    } else {
        // Local alignment: search entire M matrix
        max_val = T::zero();
        let m_matrix = &align_grid.m_matrix;

        // Pre-allocate with reasonable capacity
        max_loc.reserve(16);

        for row in 0..m_matrix.nrow {
            for col in 0..m_matrix.ncol {
                let val = m_matrix.get_score(row, col);
                if val > max_val && !T::fuzzy_equals(val, max_val) {
                    max_val = val;
                    max_loc.clear();
                    max_loc.insert((MatrixType::M, row, col));
                } else if T::fuzzy_equals(val, max_val) {
                    max_loc.insert((MatrixType::M, row, col));
                }
            }
        }
    }

    (max_val, max_loc)
}

/// Perform traceback from a specific position using parent pointers
fn traceback_from_position<T: Clone + Copy + Display + FromStr + Zero>(
    align_grid: &AlignGrid<T>,
    alignment_parameters: &AlignmentParameters<T>,
    writer: &mut BufWriter<File>,
    start_matrix: MatrixType,
    start_row: usize,
    start_col: usize,
) -> Result<(), Box<dyn Error>> {
    // Store nodes with parent indices to avoid cloning paths
    struct Node {
        pointer: Pointer,
        parent: Option<usize>,
    }

    let mut nodes = Vec::with_capacity(1024);
    let mut stack = Vec::with_capacity(256);
    let mut leaf_nodes = Vec::with_capacity(256);

    // Add initial node
    nodes.push(Node {
        pointer: (start_matrix, start_row, start_col),
        parent: None,
    });
    stack.push(0);

    // Build traceback tree using iterative DFS
    while let Some(node_idx) = stack.pop() {
        let (matrix, row, col) = nodes[node_idx].pointer;

        let pointers = match matrix {
            M => align_grid.m_matrix.get_pointers(row, col),
            Ix => align_grid.ix_matrix.get_pointers(row, col),
            Iy => align_grid.iy_matrix.get_pointers(row, col),
        };

        if pointers.is_empty() {
            leaf_nodes.push(node_idx);
            continue;
        }

        for &pointer in pointers {
            let new_node_idx = nodes.len();
            nodes.push(Node {
                pointer,
                parent: Some(node_idx),
            });
            stack.push(new_node_idx);
        }
    }

    // Process in chunks to balance memory and I/O

    let seq_a_chars = &alignment_parameters.sequences.seq_a;
    let seq_b_chars = &alignment_parameters.sequences.seq_b;

    // Pre-allocate alignment buffers
    let estimated_len = seq_a_chars.len().max(seq_b_chars.len()) + 100;

    for chunk in leaf_nodes.chunks(CHUNK_SIZE) {
        for &leaf_idx in chunk {
            // Reconstruct path by walking back through parents
            let mut path_len = 0;
            let mut current_idx = Some(leaf_idx);
            while let Some(idx) = current_idx {
                path_len += 1;
                current_idx = nodes[idx].parent;
            }

            // Build alignment strings directly
            let mut align_a = String::with_capacity(estimated_len);
            let mut align_b = String::with_capacity(estimated_len);

            current_idx = Some(leaf_idx);
            while let Some(idx) = current_idx {
                let (m, r, c) = nodes[idx].pointer;

                match m {
                    M => {
                        align_a.push(seq_a_chars[r]);
                        align_b.push(seq_b_chars[c]);
                    }
                    Ix => {
                        if c < align_grid.ix_matrix.ncol - 1 {
                            align_a.push(seq_a_chars[r]);
                            align_b.push('_');
                        }
                    }
                    Iy => {
                        if r < align_grid.iy_matrix.nrow - 1 {
                            align_a.push('_');
                            align_b.push(seq_b_chars[c]);
                        }
                    }
                }

                current_idx = nodes[idx].parent;
            }

            // Write alignment
            writer.write_all(b"\n")?;
            writer.write_all(align_a.as_bytes())?;
            writer.write_all(b"\n")?;
            writer.write_all(align_b.as_bytes())?;
            writer.write_all(b"\n")?;
        }
    }

    Ok(())
}

/// Perform traceback to generate alignments
pub fn traceback<T: Copy + FromStr + Display + Epsilon + PartialOrd + Zero>(
    align_grid: &AlignGrid<T>,
    alignment_parameters: &AlignmentParameters<T>,
    output_file: &str
) -> Result<(), Box<dyn Error>> {
    let (max_val, max_loc) = find_traceback_start(align_grid, alignment_parameters);

    let file = File::create(output_file)?;
    let mut writer = BufWriter::with_capacity(65536, file);

    writeln!(writer, "{}", max_val)?;

    for (matrix, row, col) in max_loc {
        traceback_from_position(
            align_grid,
            alignment_parameters,
            &mut writer,
            matrix,
            row,
            col
        )?;
    }

    writer.flush()?;
    Ok(())
}