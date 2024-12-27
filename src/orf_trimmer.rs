/* Algorithm
1) For each seq, identify all the start codons in reading frames 1,2, and 3
2) For all seqs determine the true start locus based on the mode of start codons and set the frame of the alignment
3) Identify the first stop codon for all seqs in this frame
4) For all seqs determine the true stop locus based on the mode of stop codons
5) Optionally output a report on the start and stop
6) Optionally perform trimming and output the resulting output fasta file
*/
use crate::fasta_manager::Fasta;
use std::collections::HashMap;

pub(super) fn trim_to_orf(inp_fasta: &Fasta) -> Result<Fasta, String> {
    let num_seqs = inp_fasta.get_num_entries();
    let starts = find_starts(&inp_fasta, num_seqs);
    let group_start = find_group_start(&starts.expect("failed to find start codons"));

    //let first_stops = find_stops(group_start.expect("failed to find group start codon"));

    dbg!(group_start);

    //println!("{:?}", starts);

    Err(String::from("still working"))
}

fn find_starts(inp_fasta: &Fasta, num_seqs: usize) -> Option<Vec<Vec<usize>>> {
    let mut starts: Vec<Vec<usize>> = vec![Vec::new(); num_seqs];

    for entry in inp_fasta {
        let mut index = 0;
        for codon in entry.get_sequence().windows(3) {
            if codon.to_ascii_uppercase() == b"ATG" || codon.to_ascii_uppercase() == b"AUG" {
                starts[entry.get_entry_num()].push(index);
            }
            index += 1;
        }
    }

    if starts.is_empty() {
        None
    } else {
        Some(starts)
    }
}

fn find_group_start(starts: &Vec<Vec<usize>>) -> Option<usize> {
    let mut start_scores: HashMap<usize, usize> = HashMap::new();
    for entry in starts {
        let mut this_score;
        for (i, start) in entry.iter().enumerate() {
            match i + 1 {
                //This scoring matrix is arbitrary and should be adjusted based on the quality of
                //results observed
                1 => this_score = 8,
                2 => this_score = 4,
                3 => this_score = 2,
                4 => this_score = 1,
                _ => this_score = 0,
            }

            if let Some(value) = start_scores.get_mut(start) {
                *value += this_score;
            } else {
                start_scores.insert(*start, this_score);
            }
        }
    }

    let mut max_value: Option<usize> = None;
    let mut max_key: Option<usize> = None;
    for (&key, &value) in &start_scores {
        if max_value.is_none() || value > max_value.unwrap() {
            max_value = Some(value);
            max_key = Some(key);
        }
    }

    max_key
}

#[allow(unused_imports)]
mod test {
    use super::*;
    use crate::fasta_manager::{open_fasta, Fasta};

    #[test]
    fn good_starts() {
        let fake_fasta_short: Fasta = open_fasta("../fake_short.fna").unwrap();
        let starts = find_starts(&fake_fasta_short, fake_fasta_short.get_num_entries());
        assert_eq!(
            starts.unwrap(),
            vec![
                vec![2, 5],
                vec![2],
                vec![2, 7],
                vec![2, 5],
                vec![2, 5],
                vec![2],
                vec![2],
                vec![],
                vec![0]
            ]
        );
    }

    #[test]
    fn no_starts() {
        let no_fasta: Fasta = Fasta::new("../fakeFile.fna");
        let starts = find_starts(&no_fasta, no_fasta.get_num_entries());
        assert_eq!(starts, None);
    }

    #[test]
    fn good_group_starts() {
        let fake_fasta_short: Fasta = open_fasta("../fake_short.fna").unwrap();
        let starts = find_starts(&fake_fasta_short, fake_fasta_short.get_num_entries());
        let group_start = find_group_start(&starts.unwrap());
        assert_eq!(group_start, Some(2));
    }

    #[test]
    fn no_group_starts() {
        let group_start = find_group_start(&vec![Vec::new()]);
        assert_eq!(group_start, None);
    }
}
