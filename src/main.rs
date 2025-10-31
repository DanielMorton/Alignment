use std::env;
use std::process;

mod alignment;
mod error;
mod io;
mod matrix;
mod parameters;
mod traceback;

use alignment::SequenceAligner;
use error::Result;

fn main() {
    if let Err(e) = run() {
        eprintln!("Error: {}", e);
        process::exit(1);
    }
}

fn run() -> Result<()> {
    let args: Vec<String> = env::args().collect();

    if args.len() != 3 {
        eprintln!("Usage: {} <input_file> <output_file>", args[0]);
        process::exit(1);
    }

    let input_file = &args[1];
    let output_file = &args[2];

    let params = parameters::AlignmentParameters::from_file(input_file)?;
    let mut aligner = SequenceAligner::new(params);
    let result = aligner.align();

    io::write_alignment_result(output_file, &result)?;

    Ok(())
}