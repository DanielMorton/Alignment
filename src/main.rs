mod alignment;
mod io;
mod models;
mod utils;

use crate::alignment::traceback;
use crate::io::parameters::AlignmentParameters;
use crate::models::AlignGrid;
use std::env;
use std::error::Error;

fn main() -> Result<(), Box<dyn Error>> {
    let args: Vec<String> = env::args().collect();

    if args.len() != 3 {
        eprintln!("Please specify an input file and an output file as args.");
        eprintln!("Usage: {} <input_file> <output_file>", args[0]);
        std::process::exit(1);
    }

    let input_file = &args[1];
    let output_file = &args[2];

    let parameters = AlignmentParameters::<f64>::load_from_file(input_file)?;
    let mut grid = AlignGrid::new(parameters.len_a(), parameters.len_b());
    let _ = grid.populate_score_matrices(&parameters)?;
    let _ = traceback(&grid, &parameters, output_file)?;
    Ok(())
}
