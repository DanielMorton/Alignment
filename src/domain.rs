#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Matrix {
    M,
    Ix,
    Iy,
}

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct Pointer {
    pub(crate) matrix: Matrix,
    pub(crate) row: usize,
    pub(crate) col: usize,
}

impl Pointer {
    pub(crate) const fn new(matrix: Matrix, row: usize, col: usize) -> Self {
        Self { matrix, row, col }
    }
}

#[derive(Debug, Clone)]
pub(crate) struct ScoreEntry {
    pub(crate) score: f64,
    pub(crate) pointers: Vec<Pointer>,
}

impl Default for ScoreEntry {
    fn default() -> Self {
        Self {
            score: 0.0,
            pointers: Vec::new(),
        }
    }
}