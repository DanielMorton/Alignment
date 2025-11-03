use ndarray::Array2;
use num_traits::Zero;
use std::fmt::Display;

/// Matrix type identifier
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum MatrixType {
    M,
    Ix,
    Iy,
}

impl MatrixType {
    pub fn as_str(&self) -> &'static str {
        match self {
            MatrixType::M => "M",
            MatrixType::Ix => "Ix",
            MatrixType::Iy => "Iy",
        }
    }
}

/// Type alias for pointer entries (matrix_type, row, col)
pub type Pointer = (MatrixType, usize, usize);

/// Score matrix used during the alignment process
pub struct ScoreMatrix<T> {
    pub matrix_type: MatrixType,
    pub nrow: usize,
    pub ncol: usize,
    pub scores: Array2<T>,
    pub pointers: Vec<Vec<Vec<Pointer>>>,
}

impl<T: Zero + Copy + Clone + Display> ScoreMatrix<T> {
    pub fn new(matrix_type: MatrixType, nrow: usize, ncol: usize) -> Self {
        let scores = Array2::zeros((nrow, ncol));
        let pointers = vec![vec![Vec::new(); ncol]; nrow];

        Self {
            matrix_type,
            nrow,
            ncol,
            scores,
            pointers,
        }
    }

    pub fn get_score(&self, row: usize, col: usize) -> T {
        self.scores[[row, col]]
    }

    pub fn set_score(&mut self, row: usize, col: usize, score: T) {
        self.scores[[row, col]] = score;
    }

    pub fn get_pointers(&self, row: usize, col: usize) -> &[Pointer] {
        &self.pointers[row][col]
    }

    pub fn set_pointers(&mut self, row: usize, col: usize, pointers: Vec<Pointer>) {
        self.pointers[row][col] = pointers;
    }

    /// Print scores for debugging
    #[allow(dead_code)]
    pub fn print_scores(&self) {
        println!("{}=", self.matrix_type.as_str());
        for r in 0..self.nrow {
            for c in 0..self.ncol {
                print!("{:.1} ", self.scores[[r, c]]);
            }
            println!();
        }
    }

    /// Print pointers for debugging
    #[allow(dead_code)]
    pub fn print_pointers(&self) {
        println!("{} Pointers=", self.matrix_type.as_str());
        for r in 0..self.nrow {
            for c in 0..self.ncol {
                print!("{:?} ", self.pointers[r][c]);
            }
            println!();
        }
    }
}
