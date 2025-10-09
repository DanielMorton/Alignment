use std::io::Write;
use std::io;

#[derive(Debug)]
pub struct AlignmentResult {
    pub(crate) score: f64,
    pub(crate) alignments: Vec<(String, String)>,
}

impl AlignmentResult {
    pub(crate) fn write_to<W: Write>(&self, mut writer: W) -> io::Result<()> {
        writeln!(writer, "{:.1}", self.score)?;

        for (seq_a, seq_b) in &self.alignments {
            writeln!(writer)?;
            writeln!(writer, "{}", seq_a)?;
            writeln!(writer, "{}", seq_b)?;
        }

        Ok(())
    }
}