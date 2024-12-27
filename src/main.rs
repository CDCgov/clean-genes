use clap::Parser;
use fasta_manager::open_fasta;
use orf_trimmer::trim_to_orf;
use process_args::Config;

mod fasta_manager;
mod math;
mod orf_trimmer;
mod process_args;

fn main() {
    let args = Config::parse();

    if args.get_module() == "TrimToORF" {
        println!("Activating module 'TrimToORF'");

        let inp_fasta =
            open_fasta(args.get_inp_fasta()).expect("Failed to open input fasta file: {}");
        let out_fasta = trim_to_orf(&inp_fasta, args.get_out_fasta());

        //dbg!(inp_fasta.get_numbered_entry(0));
        println!("trimmed fasta: {:#?}", out_fasta);
    }
    println!("{:?}", args);
}
