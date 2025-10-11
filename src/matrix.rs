const EPSILON: f64 = 1e-6;

pub fn fuzzy_equals(a: f64, b: f64) -> bool {
    (a - b).abs() < EPSILON
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum MatrixType {
    M,  // Match/mismatch matrix
    Ix, // Gap in X (horizontal)
    Iy, // Gap in Y (vertical)
}

/// A pointer to a previous cell in the dynamic programming matrix
pub type Pointer = (MatrixType, usize, usize);

/// Single entry in a score matrix
#[derive(Debug, Clone)]
pub struct ScoreEntry {
    score: f64,
    pointers: Vec<Pointer>,
}

impl ScoreEntry {
    pub fn new() -> Self {
        Self {
            score: 0.0,
            pointers: Vec::new(),
        }
    }

    pub fn with_score(score: f64) -> Self {
        Self {
            score,
            pointers: Vec::new(),
        }
    }

    pub fn score(&self) -> f64 {
        self.score
    }

    pub fn set_score(&mut self, score: f64) {
        self.score = score;
    }

    pub fn pointers(&self) -> &[Pointer] {
        &self.pointers
    }

    pub fn add_pointer(&mut self, pointer: Pointer) {
        self.pointers.push(pointer);
    }

    pub fn set_pointers(&mut self, pointers: Vec<Pointer>) {
        self.pointers = pointers;
    }
}

impl Default for ScoreEntry {
    fn default() -> Self {
        Self::new()
    }
}

/// 2D matrix of scores with traceback pointers
#[derive(Debug, Clone)]
pub struct ScoreMatrix {
    matrix_type: MatrixType,
    data: Vec<Vec<ScoreEntry>>,
    nrows: usize,
    ncols: usize,
}

impl ScoreMatrix {
    pub fn new(matrix_type: MatrixType, nrows: usize, ncols: usize) -> Self {
        let data = vec![vec![ScoreEntry::new(); ncols]; nrows];
        Self {
            matrix_type,
            data,
            nrows,
            ncols,
        }
    }

    pub fn matrix_type(&self) -> MatrixType {
        self.matrix_type
    }

    pub fn nrows(&self) -> usize {
        self.nrows
    }

    pub fn ncols(&self) -> usize {
        self.ncols
    }

    pub fn get(&self, row: usize, col: usize) -> &ScoreEntry {
        &self.data[row][col]
    }

    pub fn get_mut(&mut self, row: usize, col: usize) -> &mut ScoreEntry {
        &mut self.data[row][col]
    }

    pub fn score(&self, row: usize, col: usize) -> f64 {
        self.data[row][col].score()
    }

    pub fn set_score(&mut self, row: usize, col: usize, score: f64) {
        self.data[row][col].set_score(score);
    }

    pub fn pointers(&self, row: usize, col: usize) -> &[Pointer] {
        self.data[row][col].pointers()
    }

    pub fn set_pointers(&mut self, row: usize, col: usize, pointers: Vec<Pointer>) {
        self.data[row][col].set_pointers(pointers);
    }

    #[allow(dead_code)]
    pub fn print_scores(&self) {
        println!("{:?}=", self.matrix_type);
        for row in &self.data {
            for entry in row {
                print!("{:.1} ", entry.score());
            }
            println!();
        }
    }

    #[allow(dead_code)]
    pub fn print_pointers(&self) {
        println!("{:?} pointers:", self.matrix_type);
        for row in &self.data {
            for entry in row {
                print!("{:?} ", entry.pointers());
            }
            println!();
        }
    }
}

/// Container for all three scoring matrices used in alignment
#[derive(Debug)]
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

    pub fn get(&self, matrix_type: MatrixType) -> &ScoreMatrix {
        match matrix_type {
            MatrixType::M => &self.m,
            MatrixType::Ix => &self.ix,
            MatrixType::Iy => &self.iy,
        }
    }

    pub fn get_mut(&mut self, matrix_type: MatrixType) -> &mut ScoreMatrix {
        match matrix_type {
            MatrixType::M => &mut self.m,
            MatrixType::Ix => &mut self.ix,
            MatrixType::Iy => &mut self.iy,
        }
    }
}