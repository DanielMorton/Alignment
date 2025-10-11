use crate::error::Result;
use crate::matrix::{fuzzy_equals, MatrixSet, MatrixType};
use crate::parameters::{AlignmentParameters, AlignmentType};
use crate::traceback::{AlignmentResult, TracebackEngine};

/// Main sequence alignment engine using dynamic programming
pub struct SequenceAligner {
    params: AlignmentParameters,
    matrices: MatrixSet,
}

impl SequenceAligner {
    pub fn new(params: AlignmentParameters) -> Self {
        let nrows = params.seq_a.len();
        let ncols = params.seq_b.len();
        let matrices = MatrixSet::new(nrows, ncols);

        Self { params, matrices }
    }

    /// Execute the alignment algorithm
    pub fn align(&mut self) -> Result<AlignmentResult> {
        self.initialize_matrices();
        self.fill_matrices();
        Ok(self.traceback())
    }

    /// Initialize first row and column of matrices
    fn initialize_matrices(&mut self) {
        self.initialize_first_column();
        self.initialize_first_row();
    }

    fn initialize_first_column(&mut self) {
        let nrows = self.params.seq_a.len();

        for row in 0..nrows {
            let a = self.params.seq_a[row];
            let b = self.params.seq_b[0];
            let match_score = self.params.match_matrix.score(a, b);

            let score = match self.params.alignment_type {
                AlignmentType::Global => match_score,
                AlignmentType::Local => match_score.max(0.0),
            };

            self.matrices.m.set_score(row, 0, score);

            if row > 0 {
                self.update_ix(row, 0);
            }
        }
    }

    fn initialize_first_row(&mut self) {
        let ncols = self.params.seq_b.len();

        for col in 0..ncols {
            let a = self.params.seq_a[0];
            let b = self.params.seq_b[col];
            let match_score = self.params.match_matrix.score(a, b);

            let score = match self.params.alignment_type {
                AlignmentType::Global => match_score,
                AlignmentType::Local => match_score.max(0.0),
            };

            self.matrices.m.set_score(0, col, score);

            if col > 0 {
                self.update_iy(0, col);
            }
        }
    }

    /// Fill all matrices using dynamic programming
    fn fill_matrices(&mut self) {
        let nrows = self.params.seq_a.len();
        let ncols = self.params.seq_b.len();

        for row in 1..nrows {
            for col in 1..ncols {
                self.update_m(row, col);
                self.update_ix(row, col);
                self.update_iy(row, col);
            }
        }
    }

    /// Update M matrix (match/mismatch)
    fn update_m(&mut self, row: usize, col: usize) {
        let a = self.params.seq_a[row];
        let b = self.params.seq_b[col];
        let match_score = self.params.match_matrix.score(a, b);

        let m_prev = self.matrices.m.score(row - 1, col - 1);
        let ix_prev = self.matrices.ix.score(row - 1, col - 1);
        let iy_prev = self.matrices.iy.score(row - 1, col - 1);

        let best_prev = m_prev.max(ix_prev).max(iy_prev);
        let mut score = best_prev + match_score;

        if self.params.alignment_type == AlignmentType::Local {
            score = score.max(0.0);
        }

        let pointers = self.compute_m_pointers(row, col, score, match_score, m_prev, ix_prev, iy_prev);

        self.matrices.m.set_score(row, col, score);
        self.matrices.m.set_pointers(row, col, pointers);
    }

    fn compute_m_pointers(
        &self,
        row: usize,
        col: usize,
        score: f64,
        match_score: f64,
        m_prev: f64,
        ix_prev: f64,
        iy_prev: f64,
    ) -> Vec<(MatrixType, usize, usize)> {
        let mut pointers = Vec::new();

        match self.params.alignment_type {
            AlignmentType::Local => {
                if !fuzzy_equals(score, 0.0) {
                    if fuzzy_equals(score, m_prev + match_score) && !fuzzy_equals(m_prev, 0.0) {
                        pointers.push((MatrixType::M, row - 1, col - 1));
                    }
                    if fuzzy_equals(score, ix_prev + match_score) && !fuzzy_equals(ix_prev, 0.0) {
                        pointers.push((MatrixType::Ix, row - 1, col - 1));
                    }
                    if fuzzy_equals(score, iy_prev + match_score) && !fuzzy_equals(iy_prev, 0.0) {
                        pointers.push((MatrixType::Iy, row - 1, col - 1));
                    }
                }
            }
            AlignmentType::Global => {
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
        }

        pointers
    }

    /// Update Ix matrix (gap in X/horizontal)
    fn update_ix(&mut self, row: usize, col: usize) {
        let m_prev = self.matrices.m.score(row - 1, col);
        let ix_prev = self.matrices.ix.score(row - 1, col);

        let (dy, ey) = self.get_gap_penalties_y(col);

        let m_gap = m_prev - dy;
        let ix_gap = ix_prev - ey;

        let mut score = m_gap.max(ix_gap);
        if self.params.alignment_type == AlignmentType::Local {
            score = score.max(0.0);
        }

        let pointers = self.compute_ix_pointers(row, col, score, m_gap, ix_gap);

        self.matrices.ix.set_score(row, col, score);
        self.matrices.ix.set_pointers(row, col, pointers);
    }

    fn compute_ix_pointers(
        &self,
        row: usize,
        col: usize,
        score: f64,
        m_gap: f64,
        ix_gap: f64,
    ) -> Vec<(MatrixType, usize, usize)> {
        let mut pointers = Vec::new();

        match self.params.alignment_type {
            AlignmentType::Local => {
                if !fuzzy_equals(score, 0.0) {
                    if fuzzy_equals(score, m_gap) && !fuzzy_equals(m_gap, 0.0) {
                        pointers.push((MatrixType::M, row - 1, col));
                    }
                    if fuzzy_equals(score, ix_gap) && !fuzzy_equals(ix_gap, 0.0) {
                        pointers.push((MatrixType::Ix, row - 1, col));
                    }
                }
            }
            AlignmentType::Global => {
                if fuzzy_equals(score, m_gap) {
                    pointers.push((MatrixType::M, row - 1, col));
                }
                if fuzzy_equals(score, ix_gap) {
                    pointers.push((MatrixType::Ix, row - 1, col));
                }
            }
        }

        pointers
    }

    /// Update Iy matrix (gap in Y/vertical)
    fn update_iy(&mut self, row: usize, col: usize) {
        let m_prev = self.matrices.m.score(row, col - 1);
        let iy_prev = self.matrices.iy.score(row, col - 1);

        let (dx, ex) = self.get_gap_penalties_x(row);

        let m_gap = m_prev - dx;
        let iy_gap = iy_prev - ex;

        let mut score = m_gap.max(iy_gap);
        if self.params.alignment_type == AlignmentType::Local {
            score = score.max(0.0);
        }

        let pointers = self.compute_iy_pointers(row, col, score, m_gap, iy_gap);

        self.matrices.iy.set_score(row, col, score);
        self.matrices.iy.set_pointers(row, col, pointers);
    }

    fn compute_iy_pointers(
        &self,
        row: usize,
        col: usize,
        score: f64,
        m_gap: f64,
        iy_gap: f64,
    ) -> Vec<(MatrixType, usize, usize)> {
        let mut pointers = Vec::new();

        match self.params.alignment_type {
            AlignmentType::Local => {
                if !fuzzy_equals(score, 0.0) {
                    if fuzzy_equals(score, m_gap) && !fuzzy_equals(m_gap, 0.0) {
                        pointers.push((MatrixType::M, row, col - 1));
                    }
                    if fuzzy_equals(score, iy_gap) && !fuzzy_equals(iy_gap, 0.0) {
                        pointers.push((MatrixType::Iy, row, col - 1));
                    }
                }
            }
            AlignmentType::Global => {
                if fuzzy_equals(score, m_gap) {
                    pointers.push((MatrixType::M, row, col - 1));
                }
                if fuzzy_equals(score, iy_gap) {
                    pointers.push((MatrixType::Iy, row, col - 1));
                }
            }
        }

        pointers
    }

    /// Get gap penalties for Y direction (Ix matrix)
    /// For global alignment at boundaries, penalties may be zero
    fn get_gap_penalties_y(&self, col: usize) -> (f64, f64) {
        if self.params.alignment_type == AlignmentType::Global
            && col == self.matrices.ix.ncols() - 1
        {
            (0.0, 0.0)
        } else {
            (self.params.gap_penalties.dy, self.params.gap_penalties.ey)
        }
    }

    /// Get gap penalties for X direction (Iy matrix)
    /// For global alignment at boundaries, penalties may be zero
    fn get_gap_penalties_x(&self, row: usize) -> (f64, f64) {
        if self.params.alignment_type == AlignmentType::Global
            && row == self.matrices.iy.nrows() - 1
        {
            (0.0, 0.0)
        } else {
            (self.params.gap_penalties.dx, self.params.gap_penalties.ex)
        }
    }

    /// Perform traceback to extract alignments
    fn traceback(&self) -> AlignmentResult {
        let engine = TracebackEngine::new(&self.matrices, &self.params);
        engine.traceback()
    }
}