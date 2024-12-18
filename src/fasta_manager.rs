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

pub(crate) struct FastaEntry {
    defline: String,
    sequence: Vec<u8>,
}

impl FastaEntry {
    pub(crate) fn new(defline: String, sequence: Vec<u8>) -> Self {
        //do quality check first

        FastaEntry { defline, sequence }
    }
}

///Comes with the assumption that this is an alignment file
pub(crate) fn open_fasta(inp_fasta_name: &str) -> Result<Fasta, Box<dyn Error>> {
    let contents = fs::read_to_string(inp_fasta_name)?;

    let mut this_fasta = Fasta::new(inp_fasta_name);
    let mut last_defline = String::new();
    let mut last_seq: Vec<u8> = Vec::new();
    for line in contents.lines() {
        if line.starts_with('>') {
            if !last_seq.is_empty() {
                let this_entry = FastaEntry::new(last_defline.clone(), last_seq.clone());
                this_fasta.add(this_entry);
                last_seq = Vec::new();
            }
            last_defline = String::from(&line[1..]);
        } else {
            last_seq.extend(line[1..].as_bytes());
        }
    }
    let this_entry = FastaEntry::new(last_defline.clone(), last_seq);
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
}
