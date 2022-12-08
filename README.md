# Generator

This crate contains a library API and binary CLI to generate semi-random sequence pairs
for testing pairwise aligners.

It can generate various types of sequences using different error models. It
returns a (vector of) sequence pairs, or writes the pairs to a `.seq` file.

See `cargo run -- --help` for command line arguments and `src/lib.rs` for more details.
