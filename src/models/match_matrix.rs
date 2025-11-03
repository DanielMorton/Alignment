use num_traits::Zero;
use std::collections::HashMap;
use std::fs::File;
use std::io;
use std::io::{BufReader, Lines};
use std::str::FromStr;

/// Match matrix stores the scores of matches between characters
#[derive(Debug, Clone)]
pub struct MatchMatrix<T> {
    scores: HashMap<char, HashMap<char, T>>,
}

impl<T: Copy + FromStr + Zero> MatchMatrix<T> {
    pub fn new() -> Self {
        Self {
            scores: HashMap::new(),
        }
    }

    /// Updates or adds a score for a specified match
    fn set_score(&mut self, a: char, b: char, score: T) {
        self.scores
            .entry(a)
            .or_insert_with(HashMap::new)
            .insert(b, score);
    }

    pub fn read_match_matrix(lines: &mut Lines<BufReader<File>>) -> io::Result<Self>
    where
        <T as FromStr>::Err: std::fmt::Display,
    {
        let mut match_matrix = Self::new();
        for line in lines {
            let line = line?;
            let parts: Vec<&str> = line.split_whitespace().collect();

            if parts.len() < 5 {
                break;
            }

            // Format: "match" "(" a b score ")"
            let a = parts[2]
                .chars()
                .next()
                .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "Invalid character A"))?;
            let b = parts[3]
                .chars()
                .next()
                .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "Invalid character B"))?;
            let score = parts[4].parse().map_err(|e| {
                io::Error::new(io::ErrorKind::InvalidData, format!("Invalid score: {}", e))
            })?;

            match_matrix.set_score(a, b, score);
        }
        Ok(match_matrix)
    }

    /// Returns the score for a particular match
    pub fn get_score(&self, a: char, b: char) -> T {
        *self
            .scores
            .get(&a)
            .and_then(|m| m.get(&b))
            .unwrap_or(&T::zero())
    }
}
