use std::fs::File;
use std::io;
use std::io::{BufReader, Lines};

#[derive(Clone, Debug)]
pub struct Alphabet {
    alphabet: String,
    aliphabet_size: usize,
}

impl Alphabet {
    pub fn new(alphabet: String) -> Self {
        let size = alphabet.len();
        Self {
            alphabet,
            aliphabet_size: size,
        }
    }

    pub fn read_alphabet(lines: &mut Lines<BufReader<File>>) -> io::Result<Self> {
        let len_alphabet = lines
            .next()
            .ok_or_else(|| {
                io::Error::new(io::ErrorKind::InvalidData, "Missing alphabet A length")
            })??
            .parse::<usize>()
            .map_err(|e| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("Invalid alphabet A length: {}", e),
                )
            })?;
        let alphabet = lines
            .next()
            .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "Missing alphabet A"))??;

        // Check if len_alphabet matches the actual alphabet size
        if len_alphabet != alphabet.len() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "Alphabet length mismatch: expected {}, got {}",
                    len_alphabet,
                    alphabet.len()
                ),
            ));
        }

        Ok(Alphabet::new(alphabet))
    }
}
