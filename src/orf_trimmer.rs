use crate::fasta_manager::{remove_gaps, Fasta, FastaEntry};
use crate::math::mode_vec_usize;
use std::collections::HashMap;

/// The main functon of the TrimToORF module. Takes a Fasta object as input and returns a Fasta
/// object trimmed to what is determined to be the group start and stop codons
pub(crate) fn trim_to_orf(inp_fasta: &Fasta, out_fasta: &str) -> Result<Fasta, String> {
    let num_seqs = inp_fasta.num_entries();
    let starts = find_starts(inp_fasta, num_seqs).expect("failed to find start codons");
    let group_start = find_group_start(&starts).expect("failed to find group start codon");
    let first_stops =
        find_first_stops(inp_fasta, group_start).expect("fialed to find first stop codons");
    let group_stop = mode_vec_usize(&first_stops).expect("failed to find group stop codon");
    perform_trimming(inp_fasta, group_start, group_stop, out_fasta)
}

/// Identifies all start codons in all reading frames for a Fasta object
fn find_starts(inp_fasta: &Fasta, num_seqs: usize) -> Option<Vec<Vec<usize>>> {
    let mut starts: Vec<Vec<usize>> = vec![Vec::new(); num_seqs];

    for entry in inp_fasta {
        for (i, codon) in entry.sequence().to_ascii_uppercase().windows(3).enumerate() {
            if codon == b"ATG" || codon == b"AUG" {
                starts[entry.entry_num()].push(i);
            }
        }
    }

    if starts.is_empty() {
        None
    } else {
        Some(starts)
    }
}

/// Identifies the common start codon locus based on the location and consistency of available
/// start codons in the provided fasta file
fn find_group_start(starts: &Vec<Vec<usize>>) -> Option<usize> {
    let mut start_scores: HashMap<usize, usize> = HashMap::new();
    for entry in starts {
        let mut this_score;
        for (i, start) in entry.iter().enumerate() {
            this_score = match i + 1 {
                //This scoring matrix is arbitrary and should be adjusted based on the quality of
                //results observed
                1 => 8,
                2 => 4,
                3 => 2,
                4 => 1,
                _ => 0,
            };

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

/// Identifies the common stop codon locus. Uses the determined common start codon locus to define
/// the reading frame and then identifies the first stop codon for each sequence in that frame
fn find_first_stops(inp_fasta: &Fasta, group_start: usize) -> Option<Vec<usize>> {
    let mut first_stops: Vec<usize> = Vec::new();

    for entry in inp_fasta {
        if group_start >= entry.sequence().len() {
            return None;
        }

        if group_start >= entry.sequence().len() {
            return None;
        } else {
            for (codon_index, codon) in entry.sequence()[group_start..]
                .iter()
                .copied()
                .map(|b| b.to_ascii_uppercase())
                .enumerate()
                .filter(|(_, b)| *b != b'-')
                .array_chunks::<3>()
                .map(|a| (a[0].0, [a[0].1, a[1].1, a[2].1]))
            {
                if matches!(&codon, b"TAG" | b"TGA" | b"TAA" | b"UAG" | b"UGA" | b"UAA") {
                    first_stops.push(group_start + codon_index);
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

/// Does the actual trimming step, taking in the Fasta object, the group start and stop codons (the
/// locus at which to trim), and the name of the output file and returns a trimmed Fasta object
/// with a new name matching the name of the output file
fn perform_trimming(
    inp_fasta: &Fasta,
    start: usize,
    stop: usize,
    out_fasta_name: &str,
) -> Result<Fasta, String> {
    let mut trimmed_fasta = Fasta::new(out_fasta_name);

    for entry in inp_fasta {
        let mut trimmed_sequence: Vec<u8> = Vec::new();
        for (i, base) in entry.sequence().iter().enumerate() {
            if i >= start && i < stop + 3 {
                trimmed_sequence.push(*base);
            }
        }
        //if entry.entry_num() == 1 {}

        let trimmed_entry = FastaEntry::new(entry.defline(), trimmed_sequence, entry.entry_num());
        trimmed_fasta.add(trimmed_entry);
    }

    if trimmed_fasta.num_entries() == 0 {
        Err(String::from("failed to trim fasta"))
    } else {
        Ok(trimmed_fasta)
    }
}

#[allow(unused_imports)]
mod test {
    use super::*;
    use crate::fasta_manager::{open_fasta, Fasta};

    #[test]
    fn good_starts() {
        let fake_fasta_short: Fasta = open_fasta("fake_short.fna").unwrap();
        let starts = find_starts(&fake_fasta_short, fake_fasta_short.num_entries());
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
        let no_fasta: Fasta = Fasta::new("fakeFile.fna");
        let starts = find_starts(&no_fasta, no_fasta.num_entries());
        assert_eq!(starts, None);
    }

    #[test]
    fn good_group_starts() {
        let fake_fasta_short: Fasta = open_fasta("fake_short.fna").unwrap();
        let starts = find_starts(&fake_fasta_short, fake_fasta_short.num_entries());
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
        let fake_fasta_short: Fasta = open_fasta("fake_short.fna").unwrap();
        let starts = find_starts(&fake_fasta_short, fake_fasta_short.num_entries());
        let group_start = find_group_start(&starts.unwrap()).unwrap();
        let first_stops = find_first_stops(&fake_fasta_short, group_start);

        assert_eq!(first_stops, Some(vec![8, 5, 8, 8, 8, 8]));
    }

    #[test]
    fn bad_first_stop() {
        let fake_fasta_short: Fasta = open_fasta("fake_short.fna").unwrap();
        let group_start = 70;
        let first_stops = find_first_stops(&fake_fasta_short, group_start);

        assert_eq!(first_stops, None);
    }

    #[test]
    fn full_trim_small() {
        let fake_fasta_short: Fasta = open_fasta("fake_short.fna").unwrap();
        let trimmed_fasta = trim_to_orf(&fake_fasta_short, "./output.fasta").unwrap();
        for entry in &trimmed_fasta {
            match entry.entry_num() {
                0 => assert_eq!(entry.sequence(), b"ATGATGTAG"),
                1 => assert_eq!(entry.sequence(), b"ATGTGATAA"),
                2 => assert_eq!(entry.sequence(), b"ATG--ATGA"),
                3 => assert_eq!(entry.sequence(), b"atgatgtag"),
                4 => assert_eq!(entry.sequence(), b"atGAtGTAG"),
                5 => assert_eq!(entry.sequence(), b"ATGWKDTAG"),
                6 => assert_eq!(entry.sequence(), b"ATGKSMTAA"),
                7 => assert_eq!(entry.sequence(), b"NNNNNNNNN"),
                8 => assert_eq!(entry.sequence(), b"GNG--TTGA"),
                _ => panic!(),
            }
        }
    }
}
