# Sequence Alignment - Rust Implementation

A memory-safe, high-performance Rust implementation of sequence alignment algorithms using dynamic programming with affine gap penalties.

## Features

- **Global Alignment (Needleman-Wunsch)**: Aligns entire sequences end-to-end
- **Local Alignment (Smith-Waterman)**: Finds best local alignments within sequences
- **Affine Gap Penalties**: Uses separate gap opening and extension penalties for biological accuracy
- **Multiple Optimal Alignments**: Finds all co-optimal alignments, not just one
- **Generic Implementation**: Supports both `f64` and `i32` scoring types
- **Memory Safety**: Guaranteed by Rust's ownership system
- **Zero-Cost Abstractions**: High-level code with C-like performance

## Algorithm Details

The implementation uses three dynamic programming matrices:
- **M Matrix**: Match/mismatch scores
- **Ix Matrix**: Gaps in sequence A (insertions in B)
- **Iy Matrix**: Gaps in sequence B (insertions in A)

This approach correctly handles affine gap penalties where:
- `dx`/`dy`: Gap opening penalties
- `ex`/`ey`: Gap extension penalties

## Requirements

- **Rust 1.70+** (uses modern Rust features)
- **Cargo** (Rust's build system and package manager)

## Project Structure

```
sequence-alignment/
├── Cargo.toml                  # Project manifest
├── Cargo.lock                  # Dependency lock file
├── README.md                   # This file
└── src/
    ├── main.rs                 # Entry point
    ├── lib.rs                  # Library root (optional)
    ├── alignment/
    │   └── mod.rs             # Traceback implementation
    ├── io/
    │   └── parameters.rs      # Input parsing
    ├── models/
    │   ├── mod.rs
    │   ├── score_matrix.rs    # Score matrix structure
    │   ├── sequences.rs       # Sequence handling
    │   ├── alphabet.rs        # Alphabet handling
    │   ├── gap_penalties.rs   # Gap penalty parameters
    │   ├── match_matrix.rs    # Match/mismatch scores
    │   └── align_grid.rs      # Main alignment grid
    └── utils/
        └── mod.rs             # Utility functions
```

## Installation

### Using Cargo (Recommended)

```bash
# Clone or navigate to project directory
cd sequence-alignment

# Build in release mode (optimized)
cargo build --release

# The binary will be in target/release/alignment
./target/release/alignment input.txt output.txt

# Or build and run in one step
cargo run --release -- input.txt output.txt
```

### Development Build

```bash
# Build without optimizations (faster compilation)
cargo build

# Run with debug symbols
cargo run -- input.txt output.txt
```

## Cargo.toml

```toml
[package]
name = "alignment"
version = "1.0.0"
edition = "2021"
authors = ["Your Name <your.email@example.com>"]

[dependencies]
ndarray = "0.15"     # For 2D arrays
num-traits = "0.2"   # For numeric trait bounds

[profile.release]
opt-level = 3
lto = true           # Link-time optimization
codegen-units = 1    # Better optimization
```

## Usage

```bash
cargo run --release -- <input_file> <output_file>

# Or use the compiled binary
./target/release/alignment <input_file> <output_file>
```

### Input File Format

The input file must follow this exact format:

```
<sequence_A>
<sequence_B>
<alignment_type>
<dx> <ex> <dy> <ey>
<alphabet_A_length>
<alphabet_A>
<alphabet_B_length>
<alphabet_B>
<index_a> <index_b> <char_a> <char_b> <score>
<index_a> <index_b> <char_a> <char_b> <score>
...
```

**Parameters:**
- `alignment_type`: `0` for global, `1` for local
- `dx, ex, dy, ey`: Gap opening and extension penalties
- `alphabet_*_length`: Number of characters in the alphabet
- Each match score line has 5 fields:
    - `index_a`: 1-based index in alphabet A
    - `index_b`: 1-based index in alphabet B
    - `char_a`: Character from alphabet A
    - `char_b`: Character from alphabet B
    - `score`: Match/mismatch score

### Example Input File

```
ACGTACGT
ACGTAGCT
0
2.0 1.0 2.0 1.0
4
ACGT
4
ACGT
1 1 A A 1
1 2 A C -1
1 3 A G -1
1 4 A T -1
2 1 C A -1
2 2 C C 1
2 3 C G -1
2 4 C T -1
3 1 G A -1
3 2 G C -1
3 3 G G 1
3 4 G T -1
4 1 T A -1
4 2 T C -1
4 3 T G -1
4 4 T T 1
```

**Note:** The indices are 1-based and correspond to positions in the alphabets.

### Output File Format

```
<max_score>

<alignment_A_1>
<alignment_B_1>

<alignment_A_2>
<alignment_B_2>
...
```

The first line contains the maximum alignment score. Each subsequent pair of lines represents one optimal alignment, with gaps represented as underscores (`_`).

## Rust Features

This implementation leverages Rust's unique features:

### Safety Guarantees
- **No null pointers**: Uses `Option<T>` for optional values
- **No data races**: Compile-time thread safety checks
- **No buffer overflows**: Bounds checking on array access
- **No use-after-free**: Ownership system prevents dangling references

### Performance Features
- **Zero-cost abstractions**: High-level iterators compile to optimal code
- **Generic programming**: Monomorphization for specialized code per type
- **LLVM backend**: Same optimization infrastructure as Clang
- **No garbage collector**: Predictable performance with deterministic cleanup

### Modern Language Features
- **Pattern matching**: Exhaustive case analysis
- **Trait system**: Type-safe polymorphism
- **Error handling**: `Result<T, E>` for explicit error propagation
- **Algebraic data types**: `enum` for representing variants

## Performance Considerations

### Complexity

- **Time Complexity**: O(n × m) where n and m are sequence lengths
- **Space Complexity**: O(n × m) for storing matrices and pointers

### Optimization

The `--release` flag enables:
- **Optimization level 3**: Maximum speed optimizations
- **Link-time optimization (LTO)**: Cross-crate inlining
- **Single codegen unit**: Better global optimization
- **Dead code elimination**: Remove unused functions

### Benchmarking

```bash
# Add to Cargo.toml
[dev-dependencies]
criterion = "0.5"

# Create benches/alignment_bench.rs
# Run benchmarks
cargo bench
```

### Chunk Processing

The traceback phase processes alignments in chunks of 16,384 paths to balance memory usage and I/O efficiency.

## Architecture

### Module Organization

```rust
// src/main.rs - Entry point
mod alignment;
mod io;
mod models;
mod utils;

// src/models/mod.rs - Public API
pub mod score_matrix;
pub mod align_grid;
pub use score_matrix::{MatrixType, ScoreMatrix};
pub use align_grid::AlignGrid;
```

### Key Traits

- **`Copy + Display + FromStr + Zero`**: Required for numeric types
- **`Epsilon`**: Custom trait for fuzzy floating-point comparison
- **`Error`**: Standard error trait for error types

### Key Structures

```rust
// Generic over numeric type T
pub struct ScoreMatrix<T> {
    pub matrix_type: MatrixType,
    pub nrow: usize,
    pub ncol: usize,
    pub scores: Array2<T>,
    pub pointers: Vec<Vec<Vec<Pointer>>>,
}

pub struct AlignGrid<T> {
    pub m_matrix: ScoreMatrix<T>,
    pub ix_matrix: ScoreMatrix<T>,
    pub iy_matrix: ScoreMatrix<T>,
}
```

## Testing

### Basic Test

```bash
# Create a simple test input
cat > test_input.txt << 'EOF'
ACGT
ACCT
0
1.0 0.5 1.0 0.5
4
ACGT
4
ACGT
1 1 A A 1
1 2 A C -1
1 3 A G -1
1 4 A T -1
2 1 C A -1
2 2 C C 1
2 3 C G -1
2 4 C T -1
3 1 G A -1
3 2 G C -1
3 3 G G 1
3 4 G T -1
4 1 T A -1
4 2 T C -1
4 3 T G -1
4 4 T T 1
EOF

# Run alignment
cargo run --release -- test_input.txt test_output.txt

# View results
cat test_output.txt
```

### Unit Tests

```rust
// Add to src/models/align_grid.rs
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_new_grid() {
        let grid = AlignGrid::<f64>::new(10, 20);
        assert_eq!(grid.m_matrix.nrow, 10);
        assert_eq!(grid.m_matrix.ncol, 20);
    }
}

// Run tests
cargo test
```

### Integration Tests

```bash
# Create tests/integration_test.rs
cargo test --test integration_test
```

### RNA/DNA Alignment Example

```bash
cat > rna_dna_input.txt << 'EOF'
AUGCAUGC
ATGCATGC
1
0.7 0.1 1.5 0.5
4
AUGC
4
ATGC
1 1 A A 1
1 2 A T 0
1 3 A G 0
1 4 A C 0
2 1 U A -1
2 2 U T 1
2 3 U G 0
2 4 U C 0
3 1 G A 0
3 2 G T 0
3 3 G G 1
3 4 G C 0
4 1 C A -1
4 2 C T 0
4 3 C G 0
4 4 C C 1
EOF

cargo run --release -- rna_dna_input.txt rna_dna_output.txt
```

## Common Issues

### Compilation Errors

**"trait bounds were not satisfied"**
- Ensure generic types implement required traits
- Common fix: Add trait bounds like `T: Copy + Display + FromStr`

**"cannot find value `M` in this scope"**
- Import the enum: `use crate::models::score_matrix::MatrixType::*;`

### Runtime Errors

**"Missing sequence A"**
- Check input file format
- Ensure file ends with newline after match matrix

**"thread 'main' panicked at 'index out of bounds'"**
- Verify sequences aren't empty
- Check alphabet includes all characters in sequences
- Run with `RUST_BACKTRACE=1` for detailed error location

### Performance Issues

**Slow debug builds**
- Always use `--release` for performance testing
- Debug builds include extensive checks for safety

**High memory usage**
- For very long sequences, consider streaming output
- Monitor with: `cargo build --release && /usr/bin/time -v ./target/release/alignment`

## Development

### Linting and Formatting

```bash
# Format code
cargo fmt

# Run linter
cargo clippy

# Run linter with pedantic checks
cargo clippy -- -W clippy::pedantic
```

### Documentation

```bash
# Generate and open documentation
cargo doc --open

# Document private items too
cargo doc --document-private-items --open
```

### Profiling

```bash
# Install flamegraph
cargo install flamegraph

# Profile the application
cargo flamegraph -- input.txt output.txt

# Open flamegraph.svg in browser
```

## Extending the Code

### Adding New Scoring Schemes

Modify `MatchMatrix<T>` in `src/models/match_matrix.rs`:

```rust
impl<T: Copy + FromStr + Zero> MatchMatrix<T> {
    pub fn get_score_blosum62(&self, a: char, b: char) -> T {
        // Your custom scoring logic
    }
}
```

### Supporting Different Gap Models

Modify the `update_ix` and `update_iy` methods in `AlignGrid<T>` in `src/models/align_grid.rs`.

### Adding New Output Formats

Extend the `traceback` function in `src/alignment/mod.rs` to support JSON, XML, etc.

```rust
use serde::{Serialize, Deserialize};

#[derive(Serialize)]
struct AlignmentOutput {
    score: f64,
    alignments: Vec<(String, String)>,
}
```

## Cargo Commands Reference

```bash
# Build
cargo build              # Debug build
cargo build --release    # Release build with optimizations

# Run
cargo run -- args        # Run with arguments
cargo run --release      # Run optimized version

# Test
cargo test               # Run all tests
cargo test --release     # Run tests with optimizations
cargo test test_name     # Run specific test

# Clean
cargo clean              # Remove target directory

# Check
cargo check              # Fast compile check without codegen

# Update dependencies
cargo update             # Update to latest compatible versions

# Install
cargo install --path .   # Install to ~/.cargo/bin
```

## Dependencies

### ndarray
Provides efficient N-dimensional arrays, used for the scoring matrices.

```toml
ndarray = "0.15"
```

### num-traits
Provides numeric traits like `Zero`, `One`, etc. for generic programming.

```toml
num-traits = "0.2"
```

## Contributing

When contributing:
1. Follow Rust naming conventions (snake_case for functions/variables)
2. Run `cargo fmt` before committing
3. Ensure `cargo clippy` passes without warnings
4. Add tests for new functionality
5. Update documentation with `///` doc comments
6. Use `Result<T, E>` for fallible operations
7. Prefer iterator chains over explicit loops where appropriate

## License

MIT License

## References

- Needleman, S. B., & Wunsch, C. D. (1970). "A general method applicable to the search for similarities in the amino acid sequence of two proteins". *Journal of Molecular Biology*, 48(3), 443-453.
- Smith, T. F., & Waterman, M. S. (1981). "Identification of common molecular subsequences". *Journal of Molecular Biology*, 147(1), 195-197.
- Gotoh, O. (1982). "An improved algorithm for matching biological sequences". *Journal of Molecular Biology*, 162(3), 705-708.

## Additional Resources

- [The Rust Programming Language Book](https://doc.rust-lang.org/book/)
- [Rust by Example](https://doc.rust-lang.org/rust-by-example/)
- [Cargo Book](https://doc.rust-lang.org/cargo/)
- [ndarray Documentation](https://docs.rs/ndarray/)

## Authors

Daniel Morton

## Version History

- **1.0.0** (Current)
    - Initial Rust implementation
    - Global and local alignment support
    - Affine gap penalties
    - Multiple optimal alignment reporting
    - Memory-safe with zero unsafe code

---

For questions or issues, please open an issue on the project repository.

**Why Rust?**
- Memory safety without garbage collection
- Fearless concurrency (easily parallelizable in future)
- Zero-cost abstractions for high-level, performant code
- Excellent error messages and tooling
- Strong type system catches bugs at compile time