#![feature(iter_array_chunks)]
//test
use clap::Parser;
use fasta_manager::{open_fasta, write_fasta};
use orf_trimmer::trim_to_orf;
use process_args::Config;
use std::process;

mod fasta_manager;
mod math;
mod orf_trimmer;
mod process_args;

fn main() {
    let args = Config::parse();

    if args.module() == "TrimToORF" {
        eprintln!("Activating module 'TrimToORF'");

        let inp_fasta = match open_fasta(args.inp_fasta()) {
            Ok(success_fasta) => success_fasta,
            Err(err) => {
                eprintln!(
                    "\nFailed to open input fasta file, '{}', \nproducing the error: '{}'\n",
                    args.inp_fasta(),
                    err
                );
                process::exit(1);
            }
        };

        let out_fasta = match trim_to_orf(&inp_fasta, args.out_fasta()) {
            Ok(success_fasta) => success_fasta,
            Err(err) => {
                eprintln!("\nFailed to trim to ORF, producing the error: '{err}'\n");
                process::exit(1);
            }
        };

        write_fasta(&out_fasta);
    }
}
