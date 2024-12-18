use std::{error::Error, fs};

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
