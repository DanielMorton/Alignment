use crate::io::parameters::AlignmentParameters;
use crate::models::score_matrix::MatrixType::{Ix, Iy, M};
use crate::models::score_matrix::{MatrixType, Pointer};
use crate::models::AlignGrid;
use crate::utils::Epsilon;
use num_traits::Zero;
use std::collections::HashSet;
use std::error::Error;
use std::fmt::Display;
use std::fs::File;
use std::io::{self, Write};
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
        max_val = T::zero();
        let m_matrix = &align_grid.m_matrix;
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
    file: &mut File,
    start_matrix: MatrixType,
    start_row: usize,
    start_col: usize,
) -> Result<(), Box<dyn Error>> {
    // Store nodes with parent indices to avoid cloning paths
    struct Node {
        pointer: Pointer,
        parent: Option<usize>,
    }

    let mut nodes = Vec::new();
    let mut stack = Vec::new();
    let mut leaf_nodes = Vec::new();

    // Add initial node
    nodes.push(Node {
        pointer: (start_matrix, start_row, start_col),
        parent: None,
    });
    stack.push(0);

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


    const CHUNK_SIZE: usize = 16384;

    for chunk in leaf_nodes.chunks(CHUNK_SIZE) {
        let mut paths = Vec::new();
        for &leaf_idx in chunk {
            let mut path = Vec::new();
            let mut current_idx = Some(leaf_idx);

            while let Some(idx) = current_idx {
                path.push(nodes[idx].pointer);
                current_idx = nodes[idx].parent;
            }

            paths.push(path);
        }
        let alignments = traceback_paths(&paths, align_grid, alignment_parameters)?;
        for (align_a, align_b) in alignments {
            writeln!(file)?;
            writeln!(file, "{}", align_a)?;
            writeln!(file, "{}", align_b)?;
        }
    }

    Ok(())
}

fn traceback_paths<T: Copy + FromStr>(
    tracebacks: &[Vec<Pointer>],
    align_grid: &AlignGrid<T>,
    alignment_parameters: &AlignmentParameters<T>
) -> Result<Vec<(String, String)>, Box<dyn Error>> {
    let seq_a_chars = &alignment_parameters.sequences.seq_a;
    let seq_b_chars = &alignment_parameters.sequences.seq_b;
    let mut alignments = Vec::new();

    for traceback in tracebacks {
        let mut align_a = Vec::new();
        let mut align_b = Vec::new();

        for (m, r, c) in traceback.iter() {
            match m {
                M => {
                    align_a.push(seq_a_chars[*r]);
                    align_b.push(seq_b_chars[*c]);
                }
                Ix => {
                    if *c < align_grid.ix_matrix.ncol - 1 {
                        align_a.push(seq_a_chars[*r]);
                        align_b.push('_');
                    }
                }
                Iy => {
                    if *r < align_grid.iy_matrix.nrow - 1 {
                        align_a.push('_');
                        align_b.push(seq_b_chars[*c]);
                    }
                }
            }
        }

        alignments.push((align_a.into_iter().collect(), align_b.into_iter().collect()));
    }
    Ok(alignments)
}

/// Perform traceback to generate alignments
pub fn traceback<T: Copy + FromStr + Display + Epsilon + PartialOrd + Zero>(
    align_grid: &AlignGrid<T>,
    alignment_parameters: &AlignmentParameters<T>,
    output_file: &str
) -> Result<(), Box<dyn Error>> {
    let (max_val, max_loc) = find_traceback_start(align_grid, alignment_parameters);

    let mut file = File::create(output_file)?;
    writeln!(file, "{}", max_val)?;

    for (matrix, row, col) in max_loc {
        let _ = traceback_from_position(align_grid,
                                                 alignment_parameters,
                                                 &mut file,
                                                 matrix,
                                                 row,
                                                 col)?;
    }

    Ok(())
}

/// Write output to file
pub fn write_output<T: Display>(
    output_file: &str,
    max_val: T,
    alignments: &[(String, String)],
) -> io::Result<()> {
    let mut file = File::create(output_file)?;
    writeln!(file, "{}", max_val)?;

    for (align_a, align_b) in alignments {
        writeln!(file)?;
        writeln!(file, "{}", align_a)?;
        writeln!(file, "{}", align_b)?;
    }

    Ok(())
}
