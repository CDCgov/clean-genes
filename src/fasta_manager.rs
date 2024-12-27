use std::{error::Error, fmt, fs};

pub(crate) struct Fasta {
    filename: String,
    data: Vec<FastaEntry>,
}

impl Fasta {
    pub(crate) fn new(filename: &str) -> Self {
        Fasta {
            filename: String::from(filename),
            data: Vec::new(),
        }
    }

    pub(crate) fn add(&mut self, new_entry: FastaEntry) {
        self.data.push(new_entry);
    }

    pub(crate) fn get_num_entries(&self) -> usize {
        self.data.len()
    }

    pub(crate) fn get_numbered_entry(&self, num_entry: usize) -> &FastaEntry {
        &self.data[num_entry]
    }
}

impl fmt::Display for Fasta {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(f, "{}", self.filename)
    }
}

impl fmt::Debug for Fasta {
    fn fmt(&self, fmt: &mut fmt::Formatter<'_>) -> fmt::Result {
        fmt.debug_struct("Fasta")
            .field("filename", &self.filename)
            .field("data", &format_args!("{} sequences", &self.data.len()))
            .finish()
    }
}

impl IntoIterator for Fasta {
    type Item = FastaEntry; // The type of items we are iterating over
    type IntoIter = std::vec::IntoIter<FastaEntry>; // The iterator type

    fn into_iter(self) -> Self::IntoIter {
        self.data.into_iter()
    }
}

impl<'a> IntoIterator for &'a Fasta {
    type Item = &'a FastaEntry;
    type IntoIter = std::slice::Iter<'a, FastaEntry>;

    fn into_iter(self) -> Self::IntoIter {
        self.data.iter()
    }
}

pub(crate) struct FastaEntry {
    defline: String,
    sequence: Vec<u8>,
    entry_number: usize,
}

impl fmt::Display for FastaEntry {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(
            f,
            "name: {}, entry_number: {}",
            self.defline, self.entry_number
        )
    }
}

impl fmt::Debug for FastaEntry {
    fn fmt(&self, fmt: &mut fmt::Formatter<'_>) -> fmt::Result {
        let sequence_string = String::from_utf8(self.get_sequence().clone())
            .unwrap_or_else(|_| "Invalid UTF-8 in fasta entry".to_string());

        fmt.debug_struct("Fasta")
            .field("defline", &self.defline)
            .field("entry_number", &self.entry_number)
            .field("sequence_data", &sequence_string)
            .finish()
    }
}

impl FastaEntry {
    pub(crate) fn new(defline: String, sequence: Vec<u8>, entry_number: usize) -> Self {
        //do quality check first

        FastaEntry {
            defline,
            sequence,
            entry_number,
        }
    }

    pub(crate) fn get_defline(&self) -> String {
        self.defline.clone()
    }

    pub(crate) fn get_entry_num(&self) -> usize {
        self.entry_number
    }

    pub(crate) fn get_sequence(&self) -> &Vec<u8> {
        &self.sequence
    }
}

///Comes with the assumption that this is an alignment file
pub(crate) fn open_fasta(inp_fasta_name: &str) -> Result<Fasta, Box<dyn Error>> {
    let contents = fs::read_to_string(inp_fasta_name)?;

    let mut this_fasta = Fasta::new(inp_fasta_name);
    let mut last_defline = String::new();
    let mut last_seq: Vec<u8> = Vec::new();
    let mut entry_num = 0;
    for line in contents.lines() {
        if line.starts_with('>') {
            if !last_seq.is_empty() {
                let this_entry = FastaEntry::new(last_defline.clone(), last_seq.clone(), entry_num);
                this_fasta.add(this_entry);
                last_seq = Vec::new();
                entry_num += 1;
            }
            last_defline = String::from(&line[1..]);
        } else {
            last_seq.extend(line.as_bytes());
        }
    }
    let this_entry = FastaEntry::new(last_defline.clone(), last_seq, entry_num);
    this_fasta.add(this_entry);

    Ok(this_fasta)
}

#[allow(unused_imports)]
mod test {
    use super::*;

    #[test]
    fn test_fasta_open() {
        let fasta_name = "test_data/a_ha_h3_raw_500.fna";
        let fasta = open_fasta(fasta_name).unwrap();
        let num_entries = fasta.get_num_entries();
        assert_eq!(num_entries, 17);
    }

    #[ignore]
    #[test]
    fn print_fasta_entry() {
        let fasta_entry = FastaEntry::new(String::from("test"), b"ATGTTTCCCTGA".to_vec(), 1);

        println!("{}", fasta_entry);
        println!("{:?}", fasta_entry);
        println!("{:#?}", fasta_entry);
    }
}
