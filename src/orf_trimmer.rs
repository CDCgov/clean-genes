/* Algorithm
1) For each seq, identify all the start codons in reading frames 1,2, and 3
2) For all seqs determine the true start locus based on the mode of start codons and set the frame of the alignment
3) Identify the first stop codon for all seqs in this frame
4) For all seqs determine the true stop locus based on the mode of stop codons
5) Optionally output a report on the start and stop
6) Optionally perform trimming and output the resulting output fasta file
*/
use crate::fasta_manager::Fasta;
use crate::math::get_mode_vec_usize;
use std::collections::HashMap;

pub(crate) fn trim_to_orf(inp_fasta: &Fasta) -> Result<Fasta, String> {
    let num_seqs = inp_fasta.get_num_entries();
    let starts = find_starts(&inp_fasta, num_seqs).expect("failed to find start codons");
    let group_start = find_group_start(&starts).expect("failed to find group start codon");
    let first_stops =
        find_first_stops(&inp_fasta, group_start).expect("fialed to find first stop codons");
    let group_stop = get_mode_vec_usize(&first_stops).expect("failed to find group stop codon");

    dbg!(first_stops);

    //println!("{:?}", starts);

    Err(String::from("still working"))
}

fn find_starts(inp_fasta: &Fasta, num_seqs: usize) -> Option<Vec<Vec<usize>>> {
    let mut starts: Vec<Vec<usize>> = vec![Vec::new(); num_seqs];

    for entry in inp_fasta {
        for (i, codon) in entry
            .get_sequence()
            .to_ascii_uppercase()
            .windows(3)
            .enumerate()
        {
            if codon == b"ATG" || codon == b"AUG" {
                starts[entry.get_entry_num()].push(i);
            }
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

fn find_first_stops(inp_fasta: &Fasta, group_start: usize) -> Option<Vec<usize>> {
    let mut first_stops: Vec<usize> = Vec::new();

    for entry in inp_fasta {
        println!("\n");
        dbg!(entry.get_sequence());
        if group_start >= inp_fasta.get_num_entries() {
            return None;
        } else {
            for (i, codon) in entry.get_sequence()[group_start..]
                .to_ascii_uppercase()
                .chunks(3)
                .enumerate()
            {
                dbg!(codon);
                if codon == b"TAG"
                    || codon == b"TGA"
                    || codon == b"TAA"
                    || codon == b"UAG"
                    || codon == b"UGA"
                    || codon == b"UAA"
                {
                    first_stops.push(i * 3 + group_start);
                    break;
                }
            }
        }
    }

    if first_stops.is_empty() {
        None
    } else {
        Some(first_stops)
    }
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

    #[test]
    fn good_first_stops() {
        let fake_fasta_short: Fasta = open_fasta("../fake_short.fna").unwrap();
        let starts = find_starts(&fake_fasta_short, fake_fasta_short.get_num_entries());
        let group_start = find_group_start(&starts.unwrap()).unwrap();
        let first_stops = find_first_stops(&fake_fasta_short, group_start);

        assert_eq!(first_stops, Some(vec![8, 5, 8, 8, 8, 8, 8, 8]));
    }

    #[test]
    fn bad_first_stop() {
        let fake_fasta_short: Fasta = open_fasta("../fake_short.fna").unwrap();
        let group_start = 70;
        let first_stops = find_first_stops(&fake_fasta_short, group_start);

        assert_eq!(first_stops, None);
    }
}
