use std::collections::HashSet;

use crate::matrix::{fuzzy_equals, MatrixSet, MatrixType, Pointer};
use crate::parameters::{AlignmentParameters, AlignmentType};

/// Result of an alignment operation
#[derive(Debug)]
pub struct AlignmentResult {
    pub score: f64,
    pub alignments: Vec<(String, String)>,
}

/// Traceback handler for extracting optimal alignments
pub struct TracebackEngine<'a> {
    matrices: &'a MatrixSet,
    params: &'a AlignmentParameters,
}

impl<'a> TracebackEngine<'a> {
    pub fn new(matrices: &'a MatrixSet, params: &'a AlignmentParameters) -> Self {
        Self { matrices, params }
    }

    /// Find starting positions for traceback
    pub fn find_start_positions(&self) -> (f64, HashSet<Pointer>) {
        match self.params.alignment_type {
            AlignmentType::Global => self.find_global_start(),
            AlignmentType::Local => self.find_local_start(),
        }
    }

    fn find_global_start(&self) -> (f64, HashSet<Pointer>) {
        let row = self.matrices.m.nrows() - 1;
        let col = self.matrices.m.ncols() - 1;

        let m_score = self.matrices.m.score(row, col);
        let ix_score = self.matrices.ix.score(row, col);
        let iy_score = self.matrices.iy.score(row, col);

        let max_score = m_score.max(ix_score).max(iy_score);
        let mut starts = HashSet::new();

        if fuzzy_equals(m_score, max_score) {
            starts.insert((MatrixType::M, row, col));
        }
        if fuzzy_equals(ix_score, max_score) {
            starts.insert((MatrixType::Ix, row, col));
        }
        if fuzzy_equals(iy_score, max_score) {
            starts.insert((MatrixType::Iy, row, col));
        }

        (max_score, starts)
    }

    fn find_local_start(&self) -> (f64, HashSet<Pointer>) {
        let mut max_score = f64::NEG_INFINITY;
        let mut starts = HashSet::new();

        for row in 0..self.matrices.m.nrows() {
            for col in 0..self.matrices.m.ncols() {
                let score = self.matrices.m.score(row, col);

                if score > max_score {
                    max_score = score;
                    starts.clear();
                    starts.insert((MatrixType::M, row, col));
                } else if fuzzy_equals(score, max_score) {
                    starts.insert((MatrixType::M, row, col));
                }
            }
        }

        (max_score, starts)
    }

    /// Perform complete traceback from all optimal positions
    pub fn traceback(&self) -> AlignmentResult {
        let (score, starts) = self.find_start_positions();
        let mut all_paths = Vec::new();

        for &start in &starts {
            let paths = self.traceback_from(start);
            all_paths.extend(paths);
        }

        let alignments = all_paths
            .into_iter()
            .map(|path| self.path_to_alignment(&path))
            .collect();

        AlignmentResult { score, alignments }
    }

    /// Recursively trace back from a given position
    fn traceback_from(&self, position: Pointer) -> Vec<Vec<Pointer>> {
        let (matrix_type, row, col) = position;

        // Base case: reached origin
        if row == 0 && col == 0 {
            return vec![vec![position]];
        }

        let pointers = self.matrices.get(matrix_type).pointers(row, col);

        // No pointers means this is a starting position (local alignment)
        if pointers.is_empty() {
            return vec![vec![position]];
        }

        // Recursively trace back from all pointers
        let mut all_paths = Vec::new();
        for &pointer in pointers {
            for mut path in self.traceback_from(pointer) {
                path.push(position);
                all_paths.push(path);
            }
        }

        all_paths
    }

    /// Convert a traceback path to aligned sequences
    fn path_to_alignment(&self, path: &[Pointer]) -> (String, String) {
        let mut seq_a = Vec::new();
        let mut seq_b = Vec::new();

        for &(matrix_type, row, col) in path.iter().rev() {
            match matrix_type {
                MatrixType::M => {
                    seq_a.push(self.params.seq_a[row]);
                    seq_b.push(self.params.seq_b[col]);
                }
                MatrixType::Ix => {
                    // Gap in X direction (insert in sequence B)
                    if col < self.matrices.ix.ncols() - 1 {
                        seq_a.push(self.params.seq_a[row]);
                        seq_b.push('_');
                    }
                }
                MatrixType::Iy => {
                    // Gap in Y direction (insert in sequence A)
                    if row < self.matrices.iy.nrows() - 1 {
                        seq_a.push('_');
                        seq_b.push(self.params.seq_b[col]);
                    }
                }
            }
        }

        (
            seq_a.into_iter().rev().collect(),
            seq_b.into_iter().rev().collect(),
        )
    }
}