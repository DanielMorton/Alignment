use crate::io::parameters::AlignmentParameters;
use crate::models::score_matrix::MatrixType::{Ix, Iy, M};
use crate::models::score_matrix::{MatrixType, ScoreMatrix};
use crate::utils::{clamp_to_zero, max, Epsilon};
use num_traits::Zero;
use std::error::Error;
use std::fmt::Display;
use std::ops::Sub;
use std::str::FromStr;

/// Main alignment object
pub struct AlignGrid<T> {
    pub(crate) m_matrix: ScoreMatrix<T>,
    pub(crate) ix_matrix: ScoreMatrix<T>,
    pub(crate) iy_matrix: ScoreMatrix<T>,
}

impl<T: Copy + Display + Epsilon + FromStr + PartialEq + PartialOrd + Sub<Output = T> + Zero>
    AlignGrid<T>
{
    pub fn new(rows: usize, columns: usize) -> Self {
        Self {
            m_matrix: ScoreMatrix::new(M, rows, columns),
            ix_matrix: ScoreMatrix::new(Ix, rows, columns),
            iy_matrix: ScoreMatrix::new(Iy, rows, columns),
        }
    }

    /// Populate the score matrices
    pub fn populate_score_matrices(
        &mut self,
        alignment_parameters: &AlignmentParameters<T>,
    ) -> Result<(), Box<dyn Error>> {
        let sequences = &alignment_parameters.sequences;
        let match_matrix = &alignment_parameters.match_matrix;
        let (rows, columns) = (sequences.len_a(), sequences.len_b());

        self.m_matrix = ScoreMatrix::new(MatrixType::M, rows, columns);
        self.ix_matrix = ScoreMatrix::new(MatrixType::Ix, rows, columns);
        self.iy_matrix = ScoreMatrix::new(MatrixType::Iy, rows, columns);

        let seq_a_chars = &sequences.seq_a;
        let seq_b_chars = &sequences.seq_b;
        let is_global = alignment_parameters.global_alignment;

        // Initialize first column
        for r in 0..rows {
            let score = match_matrix.get_score(seq_a_chars[r], seq_b_chars[0]);
            let score = if !is_global {
                clamp_to_zero(score)
            } else {
                score
            };
            self.m_matrix.set_score(r, 0, score);
            if r > 0 {
                self.update_ix(&alignment_parameters, r, 0);
            }
        }

        // Initialize first row
        for c in 0..columns {
            let score = match_matrix.get_score(seq_a_chars[0], seq_b_chars[c]);
            let score = if !is_global {
                clamp_to_zero(score)
            } else {
                score
            };
            self.m_matrix.set_score(0, c, score);
            if c > 0 {
                self.update_iy(&alignment_parameters, 0, c);
            }
        }

        // Fill the rest of the matrix
        for r in 1..rows {
            for c in 1..columns {
                self.update(&alignment_parameters, r, c);
            }
        }
        Ok(())
    }

    /// Update all matrices at a given position
    fn update(&mut self, alignment_parameters: &AlignmentParameters<T>, row: usize, col: usize) {
        self.update_m(alignment_parameters, row, col);
        self.update_ix(alignment_parameters, row, col);
        self.update_iy(alignment_parameters, row, col);
    }

    /// Update M matrix at position
    fn update_m(&mut self, alignment_parameters: &AlignmentParameters<T>, row: usize, col: usize) {
        let sequences = &alignment_parameters.sequences;
        let seq_a_chars = &sequences.seq_a;
        let seq_b_chars = &sequences.seq_b;
        let score = alignment_parameters
            .match_matrix
            .get_score(seq_a_chars[row], seq_b_chars[col]);

        let m = self.m_matrix.get_score(row - 1, col - 1);
        let ix = self.ix_matrix.get_score(row - 1, col - 1);
        let iy = self.iy_matrix.get_score(row - 1, col - 1);

        let mut new_score = max(max(m, ix), iy) + score;
        if !alignment_parameters.global_alignment {
            new_score = clamp_to_zero(new_score);
        }

        let mut pointers = Vec::new();
        if !alignment_parameters.global_alignment {
            if new_score > T::epsilon() {
                if T::fuzzy_equals(new_score, m + score) && m > T::epsilon() {
                    pointers.push((M, row - 1, col - 1));
                }
                if T::fuzzy_equals(new_score, ix + score) && ix > T::epsilon() {
                    pointers.push((Ix, row - 1, col - 1));
                }
                if T::fuzzy_equals(new_score, iy + score) && iy > T::epsilon() {
                    pointers.push((Iy, row - 1, col - 1));
                }
            }
        } else {
            if T::fuzzy_equals(new_score, m + score) {
                pointers.push((M, row - 1, col - 1));
            }
            if T::fuzzy_equals(new_score, ix + score) {
                pointers.push((Ix, row - 1, col - 1));
            }
            if T::fuzzy_equals(new_score, iy + score) {
                pointers.push((Iy, row - 1, col - 1));
            }
        }

        self.m_matrix.set_score(row, col, new_score);
        self.m_matrix.set_pointers(row, col, pointers);
    }

    /// Update Ix matrix at position
    fn update_ix(&mut self, alignment_parameters: &AlignmentParameters<T>, row: usize, col: usize) {
        let mut pointers = Vec::new();
        let new_score;

        if !alignment_parameters.global_alignment {
            let m = self.m_matrix.get_score(row - 1, col) - alignment_parameters.gap_penalties.dy;
            let ix = self.ix_matrix.get_score(row - 1, col) - alignment_parameters.gap_penalties.ey;
            new_score = clamp_to_zero(max(m, ix));

            if new_score > T::epsilon() {
                if T::fuzzy_equals(new_score, m) && m > T::epsilon() {
                    pointers.push((M, row - 1, col));
                }
                if T::fuzzy_equals(new_score, ix) && ix > T::epsilon() {
                    pointers.push((Ix, row - 1, col));
                }
            }
        } else {
            let ncol = self.m_matrix.ncol;
            let dy = if col < ncol - 1 {
                alignment_parameters.gap_penalties.dy
            } else {
                T::zero()
            };
            let ey = if col < ncol - 1 {
                alignment_parameters.gap_penalties.ey
            } else {
                T::zero()
            };

            let m = self.m_matrix.get_score(row - 1, col) - dy;
            let ix = self.ix_matrix.get_score(row - 1, col) - ey;
            new_score = max(m, ix);

            if T::fuzzy_equals(new_score, m) {
                pointers.push((M, row - 1, col));
            }
            if T::fuzzy_equals(new_score, ix) {
                pointers.push((Ix, row - 1, col));
            }
        }

        self.ix_matrix.set_score(row, col, new_score);
        self.ix_matrix.set_pointers(row, col, pointers);
    }

    /// Update Iy matrix at position
    fn update_iy(&mut self, alignment_parameters: &AlignmentParameters<T>, row: usize, col: usize) {
        let mut pointers = Vec::new();
        let new_score;

        if !alignment_parameters.global_alignment {
            let m = self.m_matrix.get_score(row, col - 1) - alignment_parameters.gap_penalties.dx;
            let iy = self.iy_matrix.get_score(row, col - 1) - alignment_parameters.gap_penalties.ex;
            new_score = clamp_to_zero(max(m, iy));

            if new_score > T::epsilon() {
                if T::fuzzy_equals(new_score, m) && m > T::epsilon() {
                    pointers.push((M, row, col - 1));
                }
                if T::fuzzy_equals(new_score, iy) && iy > T::epsilon() {
                    pointers.push((Iy, row, col - 1));
                }
            }
        } else {
            let nrow = self.m_matrix.nrow;
            let dx = if row < nrow - 1 {
                alignment_parameters.gap_penalties.dx
            } else {
                T::zero()
            };
            let ex = if row < nrow - 1 {
                alignment_parameters.gap_penalties.ex
            } else {
                T::zero()
            };

            let m = self.m_matrix.get_score(row, col - 1) - dx;
            let iy = self.iy_matrix.get_score(row, col - 1) - ex;
            new_score = max(m, iy);

            if T::fuzzy_equals(new_score, m) {
                pointers.push((M, row, col - 1));
            }
            if T::fuzzy_equals(new_score, iy) {
                pointers.push((Iy, row, col - 1));
            }
        }

        self.iy_matrix.set_score(row, col, new_score);
        self.iy_matrix.set_pointers(row, col, pointers);
    }
}
