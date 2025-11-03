use std::fs::File;
use std::io;
use std::io::{BufReader, Lines};

#[derive(Clone, Debug)]
pub struct Sequences {
    pub seq_a: Vec<char>,
    pub seq_b: Vec<char>,
}

impl Sequences {
    pub fn from_string(seq_a: String, seq_b: String) -> Self {
        Self {
            seq_a: seq_a.chars().collect(),
            seq_b: seq_b.chars().collect(),
        }
    }

    pub fn load_sequences(lines: &mut Lines<BufReader<File>>) -> io::Result<Self> {
        let seq_a = lines
            .next()
            .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "Missing sequence A"))??;
        let seq_b = lines
            .next()
            .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "Missing sequence B"))??;
        Ok(Self::from_string(seq_a, seq_b))
    }

    pub fn len_a(&self) -> usize {
        self.seq_a.len()
    }

    pub fn len_b(&self) -> usize {
        self.seq_b.len()
    }
}
