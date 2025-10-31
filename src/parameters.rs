use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};

pub struct MatchMatrix {
    scores: HashMap<(char, char), f64>,
}

impl MatchMatrix {
    pub fn new() -> Self {
        Self {
            scores: HashMap::new(),
        }
    }

    pub fn insert(&mut self, a: char, b: char, score: f64) {
        self.scores.insert((a, b), score);
    }

    #[inline(always)]
    pub fn score(&self, a: char, b: char) -> f64 {
        self.scores.get(&(a, b)).copied().unwrap_or(0.0)
    }
}

#[derive(Clone, Copy)]
pub struct GapPenalties {
    pub dx: f64,
    pub ex: f64,
    pub dy: f64,
    pub ey: f64,
}

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum AlignmentType {
    Global,
    Local,
}

pub struct AlignmentParameters {
    pub seq_a: Vec<char>,
    pub seq_b: Vec<char>,
    pub alignment_type: AlignmentType,
    pub gap_penalties: GapPenalties,
    pub match_matrix: MatchMatrix,
}

impl AlignmentParameters {
    pub fn from_file(path: &str) -> std::io::Result<Self> {
        let file = File::open(path)?;
        let reader = BufReader::new(file);
        let mut lines = reader.lines();

        let seq_a: Vec<char> = lines.next().unwrap()?.trim().chars().collect();
        let seq_b: Vec<char> = lines.next().unwrap()?.trim().chars().collect();

        let alignment_type = if lines.next().unwrap()?.trim().parse::<i32>().unwrap() == 0 {
            AlignmentType::Global
        } else {
            AlignmentType::Local
        };

        let gap_vals: Vec<f64> = lines.next().unwrap()?
            .split_whitespace()
            .filter_map(|s| s.parse().ok())
            .collect();

        let gap_penalties = GapPenalties {
            dx: gap_vals[0],
            ex: gap_vals[1],
            dy: gap_vals[2],
            ey: gap_vals[3],
        };

        // Skip alphabet sections
        lines.next();
        lines.next();
        lines.next();
        lines.next();

        let mut match_matrix = MatchMatrix::new();
        for line in lines {
            let line = line?;
            let parts: Vec<&str> = line.split_whitespace().collect();
            if parts.len() >= 5 {
                let a = parts[2].chars().next().unwrap();
                let b = parts[3].chars().next().unwrap();
                let score: f64 = parts[4].parse().unwrap();
                match_matrix.insert(a, b, score);
            }
        }

        Ok(Self {
            seq_a,
            seq_b,
            alignment_type,
            gap_penalties,
            match_matrix,
        })
    }
}