use std::fs::File;
use std::io;
use std::io::{BufRead, BufReader};
use std::path::Path;
use crate::match_matrix::MatchMatrix;

#[derive(Debug, Clone)]
pub struct AlignmentParams {
    pub(crate) seq_a: Vec<char>,
    pub(crate) seq_b: Vec<char>,
    pub(crate) global: bool,
    pub(crate) dx: f64,
    pub(crate) ex: f64,
    pub(crate) dy: f64,
    pub(crate) ey: f64,
    pub(crate) match_matrix: MatchMatrix,
}

impl AlignmentParams {
    pub(crate) fn from_file<P: AsRef<Path>>(path: P) -> io::Result<Self> {
        let file = File::open(path)?;
        let mut lines = BufReader::new(file).lines();

        let seq_a = lines.next().ok_or(io::ErrorKind::InvalidData)??.chars().collect();
        let seq_b = lines.next().ok_or(io::ErrorKind::InvalidData)??.chars().collect();

        let global = lines.next().ok_or(io::ErrorKind::InvalidData)??
            .trim()
            .parse::<i32>()
            .map_err(|_| io::ErrorKind::InvalidData)? == 0;

        let gaps: Vec<f64> = lines.next().ok_or(io::ErrorKind::InvalidData)??
            .split_whitespace()
            .filter_map(|s| s.parse().ok())
            .collect();

        if gaps.len() < 4 {
            return Err(io::ErrorKind::InvalidData.into());
        }

        // Skip alphabet lines
        for _ in 0..4 {
            lines.next();
        }

        let match_matrix = lines
            .filter_map(Result::ok)
            .filter_map(|line| {
                let parts: Vec<_> = line.split_whitespace().collect();
                if parts.len() >= 5 {
                    Some((
                        parts[2].chars().next()?,
                        parts[3].chars().next()?,
                        parts[4].parse().ok()?,
                    ))
                } else {
                    None
                }
            })
            .fold(MatchMatrix::new(), |mut acc, (a, b, score)| {
                acc.set_score(a, b, score);
                acc
            });

        Ok(Self {
            seq_a,
            seq_b,
            global,
            dx: gaps[0],
            ex: gaps[1],
            dy: gaps[2],
            ey: gaps[3],
            match_matrix,
        })
    }

    #[inline]
    pub(crate) const fn dims(&self) -> (usize, usize) {
        (self.seq_a.len(), self.seq_b.len())
    }
}