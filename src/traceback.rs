use std::collections::HashSet;
use crate::matrix::{fuzzy_equals, MatrixSet, MatrixType, Pointer};
use crate::parameters::{AlignmentParameters, AlignmentType};

pub struct AlignmentResult {
    pub score: f64,
    pub alignments: Vec<(String, String)>,
}

pub struct TracebackEngine<'a> {
    matrices: &'a MatrixSet,
    params: &'a AlignmentParameters,
}

impl<'a> TracebackEngine<'a> {
    pub fn new(matrices: &'a MatrixSet, params: &'a AlignmentParameters) -> Self {
        Self { matrices, params }
    }

    pub fn traceback(&self) -> AlignmentResult {
        let (score, starts) = self.find_start_positions();

        // Limit number of traceback paths for memory efficiency
        //const MAX_PATHS: usize = 100;
        //let starts_limited: Vec<_> = starts.iter().take(MAX_PATHS).copied().collect();

        /*if starts.len() > MAX_PATHS {
            eprintln!("Warning: Found {} optimal starting positions, limiting to {}",
                      starts.len(), MAX_PATHS);
        }*/

        let mut all_paths = Vec::new();

        for (idx, &start) in starts.iter().enumerate() {
            if idx % 10 == 0 && idx > 0 {
                eprintln!("Traceback progress: {}/{}", idx, starts.len());
            }

            let mut paths = self.traceback_from(start);

            // Limit paths per starting position
            /*if paths.len() > 10 {
                paths.truncate(10);
            }*/

            all_paths.extend(paths);

            // Global limit on total paths
            /*if all_paths.len() > 1000 {
                eprintln!("Warning: Limiting total paths to 1000");
                all_paths.truncate(1000);
                break;
            }*/
        }

        eprintln!("Found {} alignment paths, converting to sequences...", all_paths.len());

        let alignments = all_paths
            .into_iter()
            .map(|path| self.path_to_alignment(&path))
            .collect();

        AlignmentResult { score, alignments }
    }

    fn find_start_positions(&self) -> (f64, HashSet<Pointer>) {
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
        let mut max_score = 0.0;
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

    fn traceback_from(&self, position: Pointer) -> Vec<Vec<Pointer>> {
        let (matrix_type, row, col) = position;

        if row == 0 && col == 0 {
            return vec![vec![position]];
        }

        let pointers = self.compute_pointers(matrix_type, row, col);

        if pointers.is_empty() {
            return vec![vec![position]];
        }

        let mut all_paths = Vec::new();

        for pointer in pointers {
            let mut paths = self.traceback_from(pointer);

            // Limit exponential explosion of paths
            /*if paths.len() > 10 {
                paths.truncate(10);
            }*/

            for mut path in paths {
                path.push(position);
                all_paths.push(path);

                // Early exit if too many paths
                /*if all_paths.len() > 100 {
                    return all_paths;
                }*/
            }
        }

        all_paths
    }

    fn compute_pointers(&self, matrix_type: MatrixType, row: usize, col: usize) -> Vec<Pointer> {
        let score = self.matrices.get(matrix_type).score(row, col);
        let mut pointers = Vec::new();

        match matrix_type {
            MatrixType::M => {
                if row == 0 || col == 0 {
                    return pointers;
                }

                let match_score = self.params.match_matrix.score(
                    self.params.seq_a[row],
                    self.params.seq_b[col]
                );

                let m_prev = self.matrices.m.score(row - 1, col - 1);
                let ix_prev = self.matrices.ix.score(row - 1, col - 1);
                let iy_prev = self.matrices.iy.score(row - 1, col - 1);

                if fuzzy_equals(score, m_prev + match_score) {
                    pointers.push((MatrixType::M, row - 1, col - 1));
                }
                if fuzzy_equals(score, ix_prev + match_score) {
                    pointers.push((MatrixType::Ix, row - 1, col - 1));
                }
                if fuzzy_equals(score, iy_prev + match_score) {
                    pointers.push((MatrixType::Iy, row - 1, col - 1));
                }
            }
            MatrixType::Ix => {
                if row == 0 {
                    return pointers;
                }

                let (dy, ey) = if self.params.alignment_type == AlignmentType::Global
                    && col == self.matrices.ix.ncols() - 1
                {
                    (0.0, 0.0)
                } else {
                    (self.params.gap_penalties.dy, self.params.gap_penalties.ey)
                };

                let m_prev = self.matrices.m.score(row - 1, col);
                let ix_prev = self.matrices.ix.score(row - 1, col);

                if fuzzy_equals(score, m_prev - dy) {
                    pointers.push((MatrixType::M, row - 1, col));
                }
                if fuzzy_equals(score, ix_prev - ey) {
                    pointers.push((MatrixType::Ix, row - 1, col));
                }
            }
            MatrixType::Iy => {
                if col == 0 {
                    return pointers;
                }

                let (dx, ex) = if self.params.alignment_type == AlignmentType::Global
                    && row == self.matrices.iy.nrows() - 1
                {
                    (0.0, 0.0)
                } else {
                    (self.params.gap_penalties.dx, self.params.gap_penalties.ex)
                };

                let m_prev = self.matrices.m.score(row, col - 1);
                let iy_prev = self.matrices.iy.score(row, col - 1);

                if fuzzy_equals(score, m_prev - dx) {
                    pointers.push((MatrixType::M, row, col - 1));
                }
                if fuzzy_equals(score, iy_prev - ex) {
                    pointers.push((MatrixType::Iy, row, col - 1));
                }
            }
        }

        pointers
    }

    fn path_to_alignment(&self, path: &[Pointer]) -> (String, String) {
        let mut seq_a = Vec::with_capacity(path.len());
        let mut seq_b = Vec::with_capacity(path.len());

        for &(matrix_type, row, col) in path.iter().rev() {
            match matrix_type {
                MatrixType::M => {
                    seq_a.push(self.params.seq_a[row]);
                    seq_b.push(self.params.seq_b[col]);
                }
                MatrixType::Ix => {
                    if col < self.matrices.ix.ncols() - 1 {
                        seq_a.push(self.params.seq_a[row]);
                        seq_b.push('_');
                    }
                }
                MatrixType::Iy => {
                    if row < self.matrices.iy.nrows() - 1 {
                        seq_a.push('_');
                        seq_b.push(self.params.seq_b[col]);
                    }
                }
            }
        }

        (seq_a.into_iter().rev().collect(), seq_b.into_iter().rev().collect())
    }
}