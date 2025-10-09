use std::cmp::Ordering;
use crate::align_result::AlignmentResult;
use crate::alignment_params::AlignmentParams;
use crate::domain::{Matrix, Pointer};
use crate::score_matrix::ScoreMatrix;
use crate::utils::{fuzzy_equals, MaxScore};

pub(crate) struct AlignmentEngine {
    params: AlignmentParams,
    m: ScoreMatrix,
    ix: ScoreMatrix,
    iy: ScoreMatrix,
}

impl AlignmentEngine {
    pub(crate) fn new(params: AlignmentParams) -> Self {
        let (nrow, ncol) = params.dims();
        Self {
            params,
            m: ScoreMatrix::new(nrow, ncol),
            ix: ScoreMatrix::new(nrow, ncol),
            iy: ScoreMatrix::new(nrow, ncol),
        }
    }

    pub(crate) fn align(mut self) -> AlignmentResult {
        self.initialize_borders();
        self.fill_matrices();
        self.traceback()
    }

    fn initialize_borders(&mut self) {
        let (nrow, ncol) = self.params.dims();

        // First column
        (0..nrow).for_each(|r| {
            let score = self.params.match_matrix.get_score(
                self.params.seq_a[r],
                self.params.seq_b[0],
            );
            let score = if self.params.global { score } else { score.max(0.0) };
            self.m.set_score(r, 0, score);

            if r > 0 {
                self.update_ix(r, 0);
            }
        });

        // First row
        (0..ncol).for_each(|c| {
            let score = self.params.match_matrix.get_score(
                self.params.seq_a[0],
                self.params.seq_b[c],
            );
            let score = if self.params.global { score } else { score.max(0.0) };
            self.m.set_score(0, c, score);

            if c > 0 {
                self.update_iy(0, c);
            }
        });
    }

    fn fill_matrices(&mut self) {
        let (nrow, ncol) = self.params.dims();

        (1..nrow).for_each(|r| {
            (1..ncol).for_each(|c| {
                self.update_m(r, c);
                self.update_ix(r, c);
                self.update_iy(r, c);
            });
        });
    }

    fn update_m(&mut self, row: usize, col: usize) {
        let match_score = self.params.match_matrix.get_score(
            self.params.seq_a[row],
            self.params.seq_b[col],
        );

        let scores = [
            self.m.get_score(row - 1, col - 1),
            self.ix.get_score(row - 1, col - 1),
            self.iy.get_score(row - 1, col - 1),
        ];

        let mut new_score = scores.max_score() + match_score;
        if !self.params.global {
            new_score = new_score.max(0.0);
        }

        let pointers = self.compute_pointers(
            new_score,
            match_score,
            &scores,
            row - 1,
            col - 1,
        );

        self.m.set_score(row, col, new_score);
        self.m.set_pointers(row, col, pointers);
    }

    fn update_ix(&mut self, row: usize, col: usize) {
        let m = self.m.get_score(row - 1, col) - self.params.dy;
        let ix = self.ix.get_score(row - 1, col) - self.params.ey;

        let mut new_score = m.max(ix);
        let mut pointers = Vec::new();

        if self.params.global {
            if fuzzy_equals(new_score, m) {
                pointers.push(Pointer::new(Matrix::M, row - 1, col));
            }
            if fuzzy_equals(new_score, ix) {
                pointers.push(Pointer::new(Matrix::Ix, row - 1, col));
            }
        } else {
            new_score = new_score.max(0.0);
            if !fuzzy_equals(new_score, 0.0) {
                if fuzzy_equals(new_score, m) && !fuzzy_equals(m, 0.0) {
                    pointers.push(Pointer::new(Matrix::M, row - 1, col));
                }
                if fuzzy_equals(new_score, ix) && !fuzzy_equals(ix, 0.0) {
                    pointers.push(Pointer::new(Matrix::Ix, row - 1, col));
                }
            }
        }

        self.ix.set_score(row, col, new_score);
        self.ix.set_pointers(row, col, pointers);
    }

    fn update_iy(&mut self, row: usize, col: usize) {
        let m = self.m.get_score(row, col - 1) - self.params.dx;
        let iy = self.iy.get_score(row, col - 1) - self.params.ex;

        let mut new_score = m.max(iy);
        let mut pointers = Vec::new();

        if self.params.global {
            if fuzzy_equals(new_score, m) {
                pointers.push(Pointer::new(Matrix::M, row, col - 1));
            }
            if fuzzy_equals(new_score, iy) {
                pointers.push(Pointer::new(Matrix::Iy, row, col - 1));
            }
        } else {
            new_score = new_score.max(0.0);
            if !fuzzy_equals(new_score, 0.0) {
                if fuzzy_equals(new_score, m) && !fuzzy_equals(m, 0.0) {
                    pointers.push(Pointer::new(Matrix::M, row, col - 1));
                }
                if fuzzy_equals(new_score, iy) && !fuzzy_equals(iy, 0.0) {
                    pointers.push(Pointer::new(Matrix::Iy, row, col - 1));
                }
            }
        }

        self.iy.set_score(row, col, new_score);
        self.iy.set_pointers(row, col, pointers);
    }

    fn compute_pointers(
        &self,
        new_score: f64,
        match_score: f64,
        prev_scores: &[f64],
        row: usize,
        col: usize,
    ) -> Vec<Pointer> {
        if !self.params.global && fuzzy_equals(new_score, 0.0) {
            return Vec::new();
        }

        let matrices = [Matrix::M, Matrix::Ix, Matrix::Iy];

        matrices
            .iter()
            .zip(prev_scores)
            .filter(|&(_, &score)| {
                if self.params.global {
                    fuzzy_equals(new_score, score + match_score)
                } else {
                    fuzzy_equals(new_score, score + match_score) && !fuzzy_equals(score, 0.0)
                }
            })
            .map(|(&matrix, _)| Pointer::new(matrix, row, col))
            .collect()
    }

    fn find_start_positions(&self) -> (f64, Vec<Pointer>) {
        if self.params.global {
            self.find_global_start()
        } else {
            self.find_local_start()
        }
    }

    fn find_global_start(&self) -> (f64, Vec<Pointer>) {
        let last_col = self.m.last_col_scores()
            .map(|(r, score)| (score, Pointer::new(Matrix::M, r, self.m.ncol - 1)));

        let last_row = self.m.last_row_scores()
            .map(|(c, score)| (score, Pointer::new(Matrix::M, self.m.nrow - 1, c)));

        last_col.chain(last_row)
            .fold((f64::NEG_INFINITY, Vec::new()), |(max_score, mut positions), (score, ptr)| {
                match score.partial_cmp(&max_score) {
                    Some(Ordering::Greater) => (score, vec![ptr]),
                    Some(Ordering::Equal) => {
                        positions.push(ptr);
                        (max_score, positions)
                    }
                    _ => (max_score, positions),
                }
            })
    }

    fn find_local_start(&self) -> (f64, Vec<Pointer>) {
        self.m.iter_scores()
            .map(|(r, c, score)| (score, Pointer::new(Matrix::M, r, c)))
            .fold((f64::NEG_INFINITY, Vec::new()), |(max_score, mut positions), (score, ptr)| {
                match score.partial_cmp(&max_score) {
                    Some(Ordering::Greater) => (score, vec![ptr]),
                    Some(Ordering::Equal) => {
                        positions.push(ptr);
                        (max_score, positions)
                    }
                    _ => (max_score, positions),
                }
            })
    }

    fn traceback(self) -> AlignmentResult {
        let (max_score, starts) = self.find_start_positions();

        let alignments = starts
            .into_iter()
            .flat_map(|ptr| self.traceback_from(ptr))
            .map(|path| self.path_to_alignment(path))
            .collect();

        AlignmentResult {
            score: max_score,
            alignments: alignments,
        }
    }

    fn traceback_from(&self, ptr: Pointer) -> Vec<Vec<Pointer>> {
        if ptr.row == 0 && ptr.col == 0 {
            return vec![vec![ptr]];
        }

        let pointers = match ptr.matrix {
            Matrix::M => self.m.pointers(ptr.row, ptr.col),
            Matrix::Ix => self.ix.pointers(ptr.row, ptr.col),
            Matrix::Iy => self.iy.pointers(ptr.row, ptr.col),
        };

        if pointers.is_empty() {
            return vec![vec![ptr]];
        }

        pointers
            .iter()
            .flat_map(|p| self.traceback_from(p.clone()))
            .map(|mut path| {
                path.push(ptr.clone());
                path
            })
            .collect()
    }

    fn path_to_alignment(&self, path: Vec<Pointer>) -> (String, String) {
        path.into_iter()
            .rev()
            .fold((String::new(), String::new()), |(mut a, mut b), ptr| {
                match ptr.matrix {
                    Matrix::M => {
                        a.push(self.params.seq_a[ptr.row]);
                        b.push(self.params.seq_b[ptr.col]);
                    }
                    Matrix::Ix => {
                        a.push(self.params.seq_a[ptr.row]);
                        b.push('_');
                    }
                    Matrix::Iy => {
                        a.push('_');
                        b.push(self.params.seq_b[ptr.col]);
                    }
                }
                (a, b)
            })
            .into()
    }
}