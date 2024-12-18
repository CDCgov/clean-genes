use clap::Parser;
use fasta_manager::open_fasta;
use process_args::Config;

mod fasta_manager;
mod orf_trimmer;
mod process_args;

fn main() {
    let args = Config::parse();

    if args.get_module() == "TrimToORF" {
        println!("Activating module 'TrimToORF'");

        open_fasta(args.get_inp_fasta());
        //trim_to_orf()
    }
    println!("{:?}", args);
}
