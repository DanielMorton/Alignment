use crate::domain::{Pointer, ScoreEntry};

#[derive(Debug, Clone)]
pub struct ScoreMatrix {
    pub(crate) nrow: usize,
    pub(crate) ncol: usize,
    data: Vec<Vec<ScoreEntry>>,
}

impl ScoreMatrix {
    pub(crate) fn new(nrow: usize, ncol: usize) -> Self {
        Self {
            nrow,
            ncol,
            data: vec![vec![ScoreEntry::default(); ncol]; nrow],
        }
    }

    #[inline]
    pub(crate) fn get_score(&self, row: usize, col: usize) -> f64 {
        self.data[row][col].score
    }

    #[inline]
    pub(crate) fn set_score(&mut self, row: usize, col: usize, score: f64) {
        self.data[row][col].score = score;
    }

    #[inline]
    pub(crate) fn pointers(&self, row: usize, col: usize) -> &[crate::domain::Pointer] {
        &self.data[row][col].pointers
    }

    pub(crate) fn set_pointers(&mut self, row: usize, col: usize, pointers: Vec<Pointer>) {
        self.data[row][col].pointers = pointers;
    }

    pub(crate) fn iter_scores(&self) -> impl Iterator<Item = (usize, usize, f64)> + '_ {
        (0..self.nrow).flat_map(move |r| {
            (0..self.ncol).map(move |c| (r, c, self.get_score(r, c)))
        })
    }

    pub(crate) fn last_row_scores(&self) -> impl Iterator<Item = (usize, f64)> + '_ {
        (0..self.ncol).map(|c| (c, self.get_score(self.nrow - 1, c)))
    }

    pub(crate) fn last_col_scores(&self) -> impl Iterator<Item = (usize, f64)> + '_ {
        (0..self.nrow).map(|r| (r, self.get_score(r, self.ncol - 1)))
    }
}