use std::fmt::Display;
use std::fs::File;
use std::io;
use std::io::{BufReader, Lines};
use std::str::FromStr;

#[derive(Debug, Clone)]
pub struct GapPenalties<T> {
    pub dx: T,
    pub ex: T,
    pub dy: T,
    pub ey: T,
}

impl<T: FromStr + Copy> GapPenalties<T> {
    pub fn new(dx: T, ex: T, dy: T, ey: T) -> Self {
        GapPenalties { dx, ex, dy, ey }
    }
}

impl<T: FromStr + Copy> GapPenalties<T>
where
    <T as FromStr>::Err: Display,
{
    pub fn read_gap_penalties(lines: &mut Lines<BufReader<File>>) -> io::Result<Self> {
        let gaps_line = lines
            .next()
            .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "Missing gap penalties"))??;
        let gaps = gaps_line
            .split_whitespace()
            .map(|s| s.parse::<T>())
            .collect::<Result<Vec<_>, _>>()
            .map_err(|e| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("Invalid gap values: {}", e),
                )
            })?;

        if gaps.len() < 4 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "Not enough gap penalty values",
            ));
        }

        Ok(GapPenalties::new(gaps[0], gaps[1], gaps[2], gaps[3]))
    }
}
