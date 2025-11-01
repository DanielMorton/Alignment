const EPSILON: f64 = 1e-6;

#[inline(always)]
pub fn fuzzy_equals(a: f64, b: f64) -> bool {
    (a - b).abs() < EPSILON
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum MatrixType {
    M,
    Ix,
    Iy,
}

pub type Pointer = (MatrixType, usize, usize);

/// Flat 1D storage for better cache locality
pub struct ScoreMatrix {
    matrix_type: MatrixType,
    data: Vec<f64>,
    nrows: usize,
    ncols: usize,
}

impl ScoreMatrix {
    pub fn new(matrix_type: MatrixType, nrows: usize, ncols: usize) -> Self {
        let data = vec![0.0; nrows * ncols];
        Self {
            matrix_type,
            data,
            nrows,
            ncols,
        }
    }

    #[inline(always)]
    pub fn matrix_type(&self) -> MatrixType {
        self.matrix_type
    }

    #[inline(always)]
    pub fn nrows(&self) -> usize {
        self.nrows
    }

    #[inline(always)]
    pub fn ncols(&self) -> usize {
        self.ncols
    }

    #[inline(always)]
    fn index(&self, row: usize, col: usize) -> usize {
        row * self.ncols + col
    }

    #[inline(always)]
    pub fn score(&self, row: usize, col: usize) -> f64 {
        self.data[self.index(row, col)]
    }

    #[inline(always)]
    pub fn set_score(&mut self, row: usize, col: usize, score: f64) {
        let idx = self.index(row, col);
        self.data[idx] = score;
    }
}

pub struct MatrixSet {
    pub m: ScoreMatrix,
    pub ix: ScoreMatrix,
    pub iy: ScoreMatrix,
}

impl MatrixSet {
    pub fn new(nrows: usize, ncols: usize) -> Self {
        Self {
            m: ScoreMatrix::new(MatrixType::M, nrows, ncols),
            ix: ScoreMatrix::new(MatrixType::Ix, nrows, ncols),
            iy: ScoreMatrix::new(MatrixType::Iy, nrows, ncols),
        }
    }

    #[inline(always)]
    pub fn get(&self, matrix_type: MatrixType) -> &ScoreMatrix {
        match matrix_type {
            MatrixType::M => &self.m,
            MatrixType::Ix => &self.ix,
            MatrixType::Iy => &self.iy,
        }
    }
}