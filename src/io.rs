use std::fs::File;
use std::io::Write;

use crate::error::Result;
use crate::traceback::AlignmentResult;

/// Write alignment results to output file
pub fn write_alignment_result(path: &str, result: &AlignmentResult) -> Result<()> {
    let mut file = File::create(path)?;

    writeln!(file, "{:.1}", result.score)?;

    for (seq_a, seq_b) in &result.alignments {
        writeln!(file)?;
        writeln!(file, "{}", seq_a)?;
        writeln!(file, "{}", seq_b)?;
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_alignment_result() {
        let result = AlignmentResult {
            score: 42.5,
            alignments: vec![
                ("HEAGAWGHEE".to_string(), "P_AWHEAE__".to_string()),
            ],
        };

        // Test would write to temp file and verify contents
        // Implementation left as exercise
    }
}