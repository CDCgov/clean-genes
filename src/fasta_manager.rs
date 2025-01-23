use std::{error::Error, fmt, fs};

///Represents a fasta file. contains a filename and a vector of FastaEntrys
pub(crate) struct Fasta {
    filename: String,
    data: Vec<FastaEntry>,
}

impl Fasta {
    /// Constructor for Fasta
    pub(crate) fn new(filename: &str) -> Self {
        Fasta {
            filename: String::from(filename),
            data: Vec::new(),
        }
    }

    /// Returns the filename of the Fasta
    pub(crate) fn filename(&self) -> &str {
        &self.filename
    }

    /// Add a FastaEntry to this Fasta
    pub(crate) fn add(&mut self, new_entry: FastaEntry) {
        self.data.push(new_entry);
    }

    /// Returns the number of FastaEntrys in this Fasta
    pub(crate) fn num_entries(&self) -> usize {
        self.data.len()
    }

    /// Returns a specific FastaEntry using its position in the data vector. Position starts with 0.
    pub(crate) fn indexed_entry(&self, num_entry: usize) -> &FastaEntry {
        &self.data[num_entry]
    }
}

/// For displaying a Fasta simply. Only shows filename
impl fmt::Display for Fasta {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(f, "{}", self.filename)
    }
}

/// For detailed dispay of a Fasta
impl fmt::Debug for Fasta {
    fn fmt(&self, fmt: &mut fmt::Formatter<'_>) -> fmt::Result {
        fmt.debug_struct("Fasta")
            .field("filename", &self.filename)
            .field("data", &format_args!("{} sequences", &self.data.len()))
            .finish()
    }
}

/// Allows iteration through the FastaEntrys in Fasta
impl IntoIterator for Fasta {
    type Item = FastaEntry; // The type of items we are iterating over
    type IntoIter = std::vec::IntoIter<FastaEntry>; // The iterator type

    fn into_iter(self) -> Self::IntoIter {
        self.data.into_iter()
    }
}

/// Allows iteration thorugh the FastaEntrys in &Fasta
impl<'a> IntoIterator for &'a Fasta {
    type Item = &'a FastaEntry;
    type IntoIter = std::slice::Iter<'a, FastaEntry>;

    fn into_iter(self) -> Self::IntoIter {
        self.data.iter()
    }
}

/// Represents a defline-sequence pair from a fasta file
#[derive(Clone)]
pub(crate) struct FastaEntry {
    defline: String,
    sequence: Vec<u8>,
    entry_number: usize,
}

/// Allows simple display for FastaEntry. Shows only the defline and the
/// position number of the entry in the Fasta. Position numbers start at 0.
impl fmt::Display for FastaEntry {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(
            f,
            "name: {}, entry_number: {}",
            self.defline, self.entry_number
        )
    }
}

/// Allows detailed display for FastaEntry. Displays sequence as a String rather
/// than a vector of u8s, which is how it is stored in clean-genes
impl fmt::Debug for FastaEntry {
    fn fmt(&self, fmt: &mut fmt::Formatter<'_>) -> fmt::Result {
        let sequence_string = String::from_utf8(self.sequence().clone())
            .unwrap_or_else(|_| "Invalid UTF-8 in fasta entry".to_string());

        fmt.debug_struct("Fasta")
            .field("defline", &self.defline)
            .field("entry_number", &self.entry_number)
            .field("sequence_data", &sequence_string)
            .finish()
    }
}

impl FastaEntry {
    /// Constructor for FastaEntry
    pub(crate) fn new(defline: String, sequence: Vec<u8>, entry_number: usize) -> Self {
        //do quality check first

        FastaEntry {
            defline,
            sequence,
            entry_number,
        }
    }

    /// Returns the defline of this FastaEntry
    pub(crate) fn defline(&self) -> String {
        self.defline.clone()
    }

    /// Returns the position of this FastaEntry in the Fasta's data vector. positions start at 0.
    pub(crate) fn entry_num(&self) -> usize {
        self.entry_number
    }

    /// Returns the sequence associated with this FastaEntry
    pub(crate) fn sequence(&self) -> &Vec<u8> {
        &self.sequence
    }

    /// Prints the data contained in a FastaEntry to stdout
    pub(crate) fn print_entry(&self) {
        println!(">{}", self.defline());
        let sequence_string = String::from_utf8(self.sequence().clone()).unwrap();
        println!("{}", sequence_string)
    }
}

/// Reads a fasta file and stores it in a Fasta object.
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
    if !last_seq.is_empty() {
        let this_entry = FastaEntry::new(last_defline.clone(), last_seq, entry_num);
        this_fasta.add(this_entry);
    }

    Ok(this_fasta)
}

/// Writes a Fasta object to a fasta file
pub(crate) fn write_fasta(fasta_obj: &Fasta) {
    for entry in fasta_obj {
        entry.print_entry();
    }
}

pub(crate) fn remove_gaps(the_vec: &[u8]) -> Vec<u8> {
    the_vec
        .iter()
        .filter(|&&byte| byte != b'-')
        .cloned()
        .collect()
}

#[allow(unused_imports)]
mod test {
    use super::*;

    const FASTA_EMPTY : &str = "test_data/empty_file.fna";
    const FASTA_NAME_1 : &str = "test_data/a_ha_h3_raw_500.fna";
    const FASTA_NAME_FAKE : &str = "test_data/fake_file.fna";

    #[test]
    fn test_empty_fasta() {
        test_fasta_file(FASTA_EMPTY, 0);
    }

    #[test]
    fn test_reg_fasta() {
        let mut fasta = test_fasta_file(FASTA_NAME_1, 17);
        test_fasta_defline(&mut fasta, 0, "MW585046{A_HA_H3}");
        test_fasta_seq(&mut fasta, 0, "-----------------------------atgaagacaacca------ttattttgatactactgacccattgggcttacagtcaaaa---cccaatcaatg---acaacaacacagccacattgtgtctaggacaccatgcagtagcaaatggaacattggtaaaaacaataagtgatgatcaaattgaggtgacaaatgctacagaattagttcagagcattccaatggggaaaatatgcaacaattcgtatagaattctagat---ggaaagaattgcacattaatagatgcaatgctaggagacccccactgtgacgcctttcagtatgagaattgggacctctttatagaaagaagcagcgccttcagcaattgcta-cccatatgacatccctaactatgcatcgctccgatccattgtagcatcctcaggaacattggaattcacagcagagggattcacatggacaggtgtcactcaaaacggaagaagcggatcctgcaaaaggggatcagccgatagtttctttagccgactgaattggctaacaaaatccggaagctcttaccccacattgaatgtgacaatgcctaacaataaaaacttcgacaagctatacatctgggggatccatcacccgagctcaactaaagagcagacaaaattgtatatccaggaatcagggcgagtaacagtctcaacaaaaagaagtcaacaaacaataatccctaacattgggtctagaccatggatcagaggtcaatcaggtaggataagcatatactggaccattgtaaaacctggagatatcctaatgataaacagtaatggcaacttagttgcaccgcggggatactttaaattgaaaacagggaaaagctctgtaatgagatcagatg---tacccataga-catttgtgtgtctgaat-gtattacaccaaatggaagcatctccaacgacaagccattccaaaatgtgaacaaagttacatatggaaaatgtcccaagtatatcagacaaaacactttaaagctggccactgggatgaggaatgtaccagaaaagcaaatcagaggaatctttggggcaatagcgggattcatcgaaaacggctgggaaggaatggttgatggatggtatgggttccgataccaaaactctgaaggaacagggcaagctgcagatctaaagagcactcaagcagccatcgaccagatcaatggaaagttaaacagagtgattgaaagaaccaatgagaaattccatcaaatagagaaggaattctcagaagtagaaggaagaattcaggacttggagaaatatgtagaagacaccaaaatagacctatggtcctacaatgcagaattgctggtggctctagaaaatcaacatacaattgacttaacagatgcagaaatgaataaattgtttgagagaactagacgcctgttaagagaaaacgcagaagacatgggaggtggatgtttcaagatttaccacaaatgtaataatgcatgcattggatcaataagaaatgggacatatgaccattacatatacagagatgaagcattaaacaaccgatttcagatcaaaggtgtagagttgaaatcaggctacaaagattggatactctggatttcattcgccatatcatgcttcttaatttgcgttgttctattgggttt------------------------------------------------------------------------------------------------------");
        test_fasta_defline(&mut fasta, 16, "KY583624{A_HA_H3}");
        test_fasta_seq(&mut fasta, 16, "-----------------------------atgaagactatca------ttgctttgagctacattctatgtctggttttcgctcaaaaaattcctggaaatg---acaatagcacggcaacgctgtgccttgggcaccatgcagtaccaaacggaacgatagtgaaaacaatcacaaatg");
    }

    #[test]
    fn test_fake_fasta() {
        let mut fasta = test_fasta_file(FASTA_NAME_FAKE, 2);
        test_fasta_seq(&mut fasta, 0, "---agcataagaaaga-aga");
        test_fasta_defline(&mut fasta, 0, " fake_test 1 a ");
        test_fasta_seq(&mut fasta, 1, "--gggcta-");
        test_fasta_defline(&mut fasta, 1, " fake_test 2 a");

    }


    fn test_fasta_file(fasta_name : &str, s : usize) -> Fasta {
        let fasta = open_fasta(fasta_name).unwrap();

        assert_eq!(fasta.num_entries(), s);
        assert_eq!(fasta.filename(), fasta_name);
        fasta
    }

    fn test_fasta_seq(fasta : &mut Fasta, i : usize, seq : &str) {
        use std::str;

        let seq_orig = fasta.indexed_entry(i).sequence();

        assert_eq!(str::from_utf8(seq_orig).unwrap(), seq);
    }

    fn test_fasta_defline(fasta : &mut Fasta, i : usize, defline : &str) {
        let defline_orig = fasta.indexed_entry(i).defline();

        assert_eq!(defline_orig, defline);
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
