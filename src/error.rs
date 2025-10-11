use std::fmt;
use std::io;

#[derive(Debug)]
pub enum AlignmentError {
    Io(io::Error),
    Parse(String),
    InvalidInput(String),
}

impl fmt::Display for AlignmentError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            AlignmentError::Io(e) => write!(f, "I/O error: {}", e),
            AlignmentError::Parse(msg) => write!(f, "Parse error: {}", msg),
            AlignmentError::InvalidInput(msg) => write!(f, "Invalid input: {}", msg),
        }
    }
}

impl std::error::Error for AlignmentError {}

impl From<io::Error> for AlignmentError {
    fn from(error: io::Error) -> Self {
        AlignmentError::Io(error)
    }
}

pub type Result<T> = std::result::Result<T, AlignmentError>;