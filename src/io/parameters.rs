use crate::models::{Alphabet, GapPenalties, MatchMatrix, Sequences};
use num_traits::Zero;
use std::fmt::Display;
use std::fs::File;
use std::io;
use std::io::{BufRead, BufReader, Lines};
use std::str::FromStr;

/// Alignment parameters loaded from input file
#[derive(Debug, Clone)]
pub struct AlignmentParameters<T: FromStr + Copy> {
    pub sequences: Sequences,
    pub global_alignment: bool,
    pub gap_penalties: GapPenalties<T>,
    pub alphabet_a: Alphabet,
    pub alphabet_b: Alphabet,
    pub match_matrix: MatchMatrix<T>,
}

impl<T: Copy + FromStr + Zero> AlignmentParameters<T>
where
    <T as FromStr>::Err: Display,
{
    pub fn new(
        sequences: Sequences,
        global_alignment: bool,
        gap_penalties: GapPenalties<T>,
        alphabet_a: Alphabet,
        alphabet_b: Alphabet,
        match_matrix: MatchMatrix<T>,
    ) -> Self {
        Self {
            sequences,
            global_alignment,
            gap_penalties,
            alphabet_a,
            alphabet_b,
            match_matrix,
        }
    }

    fn read_alignment_type(lines: &mut Lines<BufReader<File>>) -> io::Result<bool> {
        let alignment_type: i32 = lines
            .next()
            .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "Missing alignment type"))??
            .parse()
            .map_err(|e| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("Invalid alignment type: {}", e),
                )
            })?;
        Ok(alignment_type == 0)
    }

    pub fn load_from_file(input_file: &str) -> io::Result<Self> {
        let file = File::open(input_file)?;
        let reader = BufReader::new(file);
        let mut lines = reader.lines();

        // Read sequences
        let sequences = Sequences::load_sequences(&mut lines)?;

        // Read alignment type
        let global_alignment = Self::read_alignment_type(&mut lines)?;

        // Read gap penalties
        let gaps = GapPenalties::<T>::read_gap_penalties(&mut lines)?;

        // Read alphabet A
        let alphabet_a = Alphabet::read_alphabet(&mut lines)?;
        // Read alphabet B
        let alphabet_b = Alphabet::read_alphabet(&mut lines)?;

        // Read match scores
        let match_matrix = MatchMatrix::<T>::read_match_matrix(&mut lines)?;

        Ok(Self::new(
            sequences,
            global_alignment,
            gaps,
            alphabet_a,
            alphabet_b,
            match_matrix,
        ))
    }

    pub fn len_a(&self) -> usize {
        self.sequences.len_a()
    }

    pub fn len_b(&self) -> usize {
        self.sequences.len_b()
    }
}
