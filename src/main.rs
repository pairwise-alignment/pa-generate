use clap::Parser;
use pa_generate::GenerateArgs;
use std::{
    io::{BufWriter, Write},
    path::PathBuf,
};

#[derive(Parser)]
#[command(author, version, about)]
struct Cli {
    /// Location of the output `.seq` file.
    output: PathBuf,

    #[clap(flatten)]
    generate_args: GenerateArgs,
}

fn main() {
    let args = Cli::parse();

    assert_eq!(
        args.output.extension().unwrap_or_default(),
        "seq",
        "Output file must have .seq extension!"
    );

    let mut f = BufWriter::new(std::fs::File::create(args.output).unwrap());
    for (a, b) in args.generate_args.generate() {
        write!(f, ">").unwrap();
        f.write_all(&a).unwrap();
        writeln!(f).unwrap();
        write!(f, "<").unwrap();
        f.write_all(&b).unwrap();
        writeln!(f).unwrap();
    }
}
