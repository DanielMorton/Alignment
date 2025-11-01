use crate::matrix::MatrixSet;
use crate::parameters::{AlignmentParameters, AlignmentType};
use crate::traceback::{AlignmentResult, TracebackEngine};

pub struct SequenceAligner {
    params: AlignmentParameters,
    matrices: MatrixSet,
}

impl SequenceAligner {
    pub fn new(params: AlignmentParameters) -> Self {
        let nrows = params.seq_a.len();
        let ncols = params.seq_b.len();

        // Check memory requirements
        let total_cells = nrows * ncols * 3; // 3 matrices
        let bytes_per_cell = 8; // f64
        let total_mb = (total_cells * bytes_per_cell) / (1024 * 1024);

        eprintln!("Allocating matrices: {}x{} cells ({} MB)", nrows, ncols, total_mb);

        let matrices = MatrixSet::new(nrows, ncols);

        Self { params, matrices }
    }

    pub fn align(&mut self) -> AlignmentResult {
        self.initialize_matrices();
        self.fill_matrices();
        self.traceback()
    }

    fn initialize_matrices(&mut self) {
        let is_local = self.params.alignment_type == AlignmentType::Local;
        let nrows = self.params.seq_a.len();
        let ncols = self.params.seq_b.len();

        // Initialize first column
        for row in 0..nrows {
            let match_score = self.params.match_matrix.score(
                self.params.seq_a[row],
                self.params.seq_b[0]
            );
            let score = if is_local { match_score.max(0.0) } else { match_score };
            self.matrices.m.set_score(row, 0, score);

            if row > 0 {
                self.update_ix(row, 0);
            }
        }

        // Initialize first row
        for col in 0..ncols {
            let match_score = self.params.match_matrix.score(
                self.params.seq_a[0],
                self.params.seq_b[col]
            );
            let score = if is_local { match_score.max(0.0) } else { match_score };
            self.matrices.m.set_score(0, col, score);

            if col > 0 {
                self.update_iy(0, col);
            }
        }
    }

    fn fill_matrices(&mut self) {
        let nrows = self.params.seq_a.len();
        let ncols = self.params.seq_b.len();

        for row in 1..nrows {
            for col in 1..ncols {
                self.update_m(row, col);
                self.update_ix(row, col);
                self.update_iy(row, col);
            }

            // Progress indicator for large alignments
            if row % 1000 == 0 {
                eprintln!("Processing row {}/{}", row, nrows);
            }
        }
    }

    #[inline]
    fn update_m(&mut self, row: usize, col: usize) {
        let match_score = self.params.match_matrix.score(
            self.params.seq_a[row],
            self.params.seq_b[col]
        );

        let m_prev = self.matrices.m.score(row - 1, col - 1);
        let ix_prev = self.matrices.ix.score(row - 1, col - 1);
        let iy_prev = self.matrices.iy.score(row - 1, col - 1);

        let best_prev = m_prev.max(ix_prev).max(iy_prev);
        let mut score = best_prev + match_score;

        if self.params.alignment_type == AlignmentType::Local && score < 0.0 {
            score = 0.0;
        }

        self.matrices.m.set_score(row, col, score);
    }

    #[inline]
    fn update_ix(&mut self, row: usize, col: usize) {
        let (dy, ey) = if self.params.alignment_type == AlignmentType::Global
            && col == self.matrices.ix.ncols() - 1
        {
            (0.0, 0.0)
        } else {
            (self.params.gap_penalties.dy, self.params.gap_penalties.ey)
        };

        let m_prev = self.matrices.m.score(row - 1, col);
        let ix_prev = self.matrices.ix.score(row - 1, col);

        let mut score = (m_prev - dy).max(ix_prev - ey);
        if self.params.alignment_type == AlignmentType::Local && score < 0.0 {
            score = 0.0;
        }

        self.matrices.ix.set_score(row, col, score);
    }

    #[inline]
    fn update_iy(&mut self, row: usize, col: usize) {
        let (dx, ex) = if self.params.alignment_type == AlignmentType::Global
            && row == self.matrices.iy.nrows() - 1
        {
            (0.0, 0.0)
        } else {
            (self.params.gap_penalties.dx, self.params.gap_penalties.ex)
        };

        let m_prev = self.matrices.m.score(row, col - 1);
        let iy_prev = self.matrices.iy.score(row, col - 1);

        let mut score = (m_prev - dx).max(iy_prev - ex);
        if self.params.alignment_type == AlignmentType::Local && score < 0.0 {
            score = 0.0;
        }

        self.matrices.iy.set_score(row, col, score);
    }

    fn traceback(&self) -> AlignmentResult {
        eprintln!("Starting traceback...");
        let engine = TracebackEngine::new(&self.matrices, &self.params);
        engine.traceback()
    }
}