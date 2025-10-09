use std::env;
use std::fs::File;
use std::io;
use crate::alignment_engine::AlignmentEngine;
use crate::alignment_params::AlignmentParams;

mod match_matrix;
mod domain;
mod utils;
mod score_matrix;
mod alignment_params;
mod alignment_engine;
mod align_result;

fn main() -> io::Result<()> {
    let args: Vec<_> = env::args().collect();

    if args.len() != 3 {
        eprintln!("Usage: {} <input_file> <output_file>", args[0]);
        std::process::exit(1);
    }

    let params = AlignmentParams::from_file(&args[1])?;
    let result = AlignmentEngine::new(params).align();
    let output = File::create(&args[2])?;
    result.write_to(output)?;

    Ok(())
}