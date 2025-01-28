use crate::fasta_manager::{Fasta, FastaEntry};
use crate::math::mode_vec_usize;
use std::collections::HashMap;
use std::fmt;

#[derive(Debug)]
pub(crate) enum OrfTrimError {
    NoStartCodons,
    NoGroupStart,
    NoStopCodons(usize),
    TrimFailed,
}

impl fmt::Display for OrfTrimError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            OrfTrimError::NoStartCodons => write!(f, "Failed to find start codons in input alignment"),
            OrfTrimError::NoGroupStart => write!(f, "Failed to find a group start codon"),
            OrfTrimError::NoStopCodons(pos) => write!(
                f,
                "Failed to find any stop codons in the frame of the group start codon at locus {pos}",
            
            ),
            OrfTrimError::TrimFailed => write!(f, "Failed to trim fasta"),
        }
    }
}

impl std::error::Error for OrfTrimError {}

/// The main functon of the TrimToORF module. Takes a Fasta object as input and
/// returns a Fasta object trimmed to what is determined to be the group start
/// and stop codons
pub(crate) fn trim_to_orf(inp_fasta: &Fasta, out_fasta: &str) -> Result<Fasta, OrfTrimError> {
    let num_seqs = inp_fasta.num_entries();
    let starts = find_starts(inp_fasta, num_seqs)?;
    let group_start = find_group_start(&starts)?;
    let first_stops = find_first_stops(inp_fasta, group_start)?;
    let group_stop = mode_vec_usize(&first_stops).map_err(|_| OrfTrimError::TrimFailed)?;
    perform_trimming(inp_fasta, group_start, group_stop, out_fasta)
}

/// Identifies all start codons in all reading frames for a Fasta object
fn find_starts(inp_fasta: &Fasta, num_seqs: usize) -> Result<Vec<Vec<usize>>, OrfTrimError> {
    let mut starts: Vec<Vec<usize>> = vec![Vec::new(); num_seqs];

    for entry in inp_fasta {
        for (i, codon) in entry.sequence().to_ascii_uppercase().windows(3).enumerate() {
            if codon == b"ATG" || codon == b"AUG" {
                starts[entry.entry_num()].push(i);
            }
        }
    }

    if starts.is_empty() {
        Err(OrfTrimError::NoStartCodons)
    } else {
        Ok(starts)
    }
}

/// Identifies the common start codon locus based on the location and
/// consistency of available start codons in the provided fasta file.
fn find_group_start(starts: &Vec<Vec<usize>>) -> Result<usize, OrfTrimError> {
    let mut start_scores: HashMap<usize, usize> = HashMap::new();
    for entry in starts {
        let mut this_score;
        for (i, start) in entry.iter().enumerate() {
            this_score = match i + 1 {
                //This scoring matrix is arbitrary and should be adjusted based
                //on the quality of results observed
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

    let mut max_value = usize::MIN;
    let mut max_key = None;
    for (&key, &value) in &start_scores {
        if value > max_value {
            max_value = value;
            max_key = Some(key);
        }
    }

    match max_key {
        Some(locus) => Ok(locus),
        None => Err(OrfTrimError::NoGroupStart),
    }
}

/// Identifies the common stop codon locus. Uses the determined common start
/// codon locus to define the reading frame and then identifies the first stop
/// codon for each sequence in that frame
fn find_first_stops(inp_fasta: &Fasta, group_start: usize) -> Result<Vec<usize>, OrfTrimError> {
    let mut first_stops: Vec<usize> = Vec::new();

    for entry in inp_fasta {
        //if the group start codon is past the length of this sequence, move to
        //the next sequence
        if group_start < entry.sequence().len() {
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
        Err(OrfTrimError::NoStopCodons(group_start + 1))
    } else {
        Ok(first_stops)
    }
}

/// Does the actual trimming step, taking in the Fasta object, the group start
/// and stop codons (the locus at which to trim), and the name of the output
/// file and returns a trimmed Fasta object with a new name matching the name
/// of the output file
fn perform_trimming(
    inp_fasta: &Fasta,
    start: usize,
    stop: usize,
    out_fasta_name: &str,
) -> Result<Fasta, OrfTrimError> {
    let mut trimmed_fasta = Fasta::new(out_fasta_name);

    for entry in inp_fasta {
        let mut trimmed_sequence: Vec<u8> = Vec::new();
        for (i, base) in entry.sequence().iter().enumerate() {
            if i >= start && i < stop + 3 {
                trimmed_sequence.push(*base);
            }
        }

        let trimmed_entry = FastaEntry::new(entry.defline(), trimmed_sequence, entry.entry_num());
        trimmed_fasta.add(trimmed_entry);
    }

    if trimmed_fasta.num_entries() == 0 {
        Err(OrfTrimError::TrimFailed)
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
            Vec::from([
                Vec::from([2, 5]),
                Vec::from([2]),
                Vec::from([2, 7]),
                Vec::from([2, 5]),
                Vec::from([2, 5]),
                Vec::from([2]),
                Vec::from([2]),
                Vec::from([]),
                Vec::from([0])
            ])
        );
    }

    #[test]
    fn no_starts() {
        let no_fasta: Fasta = Fasta::new("fakeFile.fna");
        let starts = find_starts(&no_fasta, no_fasta.num_entries());
        assert_eq!(
            starts.unwrap_err().to_string(),
            "Failed to find start codons in input alignment"
        );
    }

    #[test]
    fn good_group_starts() {
        let fake_fasta_short: Fasta = open_fasta("fake_short.fna").unwrap();
        let starts = find_starts(&fake_fasta_short, fake_fasta_short.num_entries());
        let group_start = find_group_start(&starts.unwrap());
        assert_eq!(group_start.unwrap(), 2);
    }

    #[test]
    fn no_group_starts() {
        let group_start = find_group_start(&Vec::from([Vec::new()]));
        assert_eq!(
            group_start.unwrap_err().to_string(),
            "Failed to find a group start codon"
        );
    }

    #[test]
    fn good_first_stops() {
        let fake_fasta_short: Fasta = open_fasta("fake_short.fna").unwrap();
        let starts = find_starts(&fake_fasta_short, fake_fasta_short.num_entries());
        let group_start = find_group_start(&starts.unwrap()).unwrap();
        let first_stops = find_first_stops(&fake_fasta_short, group_start);

        assert_eq!(first_stops.unwrap(), Vec::from([8, 5, 8, 8, 8, 8]));
    }

    #[test]
    fn bad_first_stop() {
        let fake_fasta_short: Fasta = open_fasta("fake_short.fna").unwrap();
        let group_start = 70;
        let first_stops = find_first_stops(&fake_fasta_short, group_start);

        assert_eq!(
            first_stops.unwrap_err().to_string(),
            "Failed to find any stop codons in the frame of the group start codon at locus 71"
        );
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
