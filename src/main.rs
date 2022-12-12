use clap::Parser;
use pa_generate::GenerateArgs;
use std::path::PathBuf;

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
    args.generate_args.generate_file(&args.output);
}
