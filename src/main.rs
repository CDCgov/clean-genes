#![feature(iter_array_chunks)]

use clap::Parser;
use fasta_manager::{open_fasta, write_fasta};
use orf_trimmer::trim_to_orf;
use process_args::Config;

mod fasta_manager;
mod math;
mod orf_trimmer;
mod process_args;

fn main() {
    let args = Config::parse();

    if args.module() == "TrimToORF" {
        eprintln!("Activating module 'TrimToORF'");

        let inp_fasta = open_fasta(args.inp_fasta()).expect("Failed to open input fasta file: {}");
        let out_fasta = trim_to_orf(&inp_fasta, args.out_fasta());
        match &out_fasta {
            Ok(success_fasta) => write_fasta(success_fasta),
            Err(_e) => println!("failed to trim to ORF"),
        }
    }
}
