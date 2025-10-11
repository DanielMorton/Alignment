use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};

use crate::error::{AlignmentError, Result};

/// Scoring matrix for character matches
#[derive(Debug, Clone)]
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

    pub fn score(&self, a: char, b: char) -> f64 {
        self.scores.get(&(a, b)).copied().unwrap_or(0.0)
    }
}

/// Gap penalty model
#[derive(Debug, Clone, Copy)]
pub struct GapPenalties {
    pub dx: f64, // Gap open penalty for sequence X (horizontal)
    pub ex: f64, // Gap extension penalty for sequence X
    pub dy: f64, // Gap open penalty for sequence Y (vertical)
    pub ey: f64, // Gap extension penalty for sequence Y
}

/// Alignment type
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AlignmentType {
    Global,
    Local,
}

/// Complete alignment parameters
#[derive(Debug, Clone)]
pub struct AlignmentParameters {
    pub seq_a: Vec<char>,
    pub seq_b: Vec<char>,
    pub alignment_type: AlignmentType,
    pub gap_penalties: GapPenalties,
    pub match_matrix: MatchMatrix,
}

impl AlignmentParameters {
    pub fn from_file(path: &str) -> Result<Self> {
        let file = File::open(path)?;
        let reader = BufReader::new(file);
        let mut lines = reader.lines();

        let seq_a = Self::read_sequence(&mut lines, "sequence A")?;
        let seq_b = Self::read_sequence(&mut lines, "sequence B")?;
        let alignment_type = Self::read_alignment_type(&mut lines)?;
        let gap_penalties = Self::read_gap_penalties(&mut lines)?;

        // Skip alphabet information (legacy format)
        Self::skip_alphabet_section(&mut lines)?;
        Self::skip_alphabet_section(&mut lines)?;

        let match_matrix = Self::read_match_matrix(&mut lines)?;

        Ok(Self {
            seq_a,
            seq_b,
            alignment_type,
            gap_penalties,
            match_matrix,
        })
    }

    fn read_sequence(
        lines: &mut impl Iterator<Item = std::io::Result<String>>,
        name: &str,
    ) -> Result<Vec<char>> {
        let seq = lines
            .next()
            .ok_or_else(|| AlignmentError::InvalidInput(format!("Missing {}", name)))??
            .trim()
            .chars()
            .collect::<Vec<_>>();
        Ok(seq)
    }

    fn read_alignment_type(
        lines: &mut impl Iterator<Item = std::io::Result<String>>,
    ) -> Result<AlignmentType> {
        let line = lines
            .next()
            .ok_or_else(|| AlignmentError::InvalidInput("Missing alignment type".into()))??;

        let value: i32 = line
            .trim()
            .parse()
            .map_err(|_| AlignmentError::Parse("Invalid alignment type".into()))?;

        Ok(if value == 0 {
            AlignmentType::Global
        } else {
            AlignmentType::Local
        })
    }

    fn read_gap_penalties(
        lines: &mut impl Iterator<Item = std::io::Result<String>>,
    ) -> Result<GapPenalties> {
        let line = lines
            .next()
            .ok_or_else(|| AlignmentError::InvalidInput("Missing gap penalties".into()))??;

        let values: Vec<f64> = line
            .split_whitespace()
            .filter_map(|s| s.parse().ok())
            .collect();

        if values.len() < 4 {
            return Err(AlignmentError::Parse(
                "Expected 4 gap penalty values".into(),
            ));
        }

        Ok(GapPenalties {
            dx: values[0],
            ex: values[1],
            dy: values[2],
            ey: values[3],
        })
    }

    fn skip_alphabet_section(
        lines: &mut impl Iterator<Item = std::io::Result<String>>,
    ) -> Result<()> {
        lines
            .next()
            .ok_or_else(|| AlignmentError::InvalidInput("Missing alphabet length".into()))??;
        lines
            .next()
            .ok_or_else(|| AlignmentError::InvalidInput("Missing alphabet".into()))??;
        Ok(())
    }

    fn read_match_matrix(
        lines: &mut impl Iterator<Item = std::io::Result<String>>,
    ) -> Result<MatchMatrix> {
        let mut matrix = MatchMatrix::new();

        for line in lines {
            let line = line?;
            let parts: Vec<&str> = line.split_whitespace().collect();

            if parts.len() >= 5 {
                if let (Some(&a_str), Some(&b_str), Some(&score_str)) =
                    (parts.get(2), parts.get(3), parts.get(4))
                {
                    let a = a_str
                        .chars()
                        .next()
                        .ok_or_else(|| AlignmentError::Parse("Invalid character A".into()))?;
                    let b = b_str
                        .chars()
                        .next()
                        .ok_or_else(|| AlignmentError::Parse("Invalid character B".into()))?;
                    let score: f64 = score_str
                        .parse()
                        .map_err(|_| AlignmentError::Parse("Invalid score".into()))?;

                    matrix.insert(a, b, score);
                }
            }
        }

        Ok(matrix)
    }
}