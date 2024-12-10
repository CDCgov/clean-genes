use clap::Parser;
use process_args::Config;

mod process_args;

fn main() {
    let args = Config::parse();

    if args.get_module() == "TrimToORF" {
        println!("Activating module 'TrimToORF'");
    }
    println!("{:?}", args);
}
