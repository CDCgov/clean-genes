use std::collections::HashMap;
use std::error::Error;

/// Calculates the mathematical mode of a vector of usizes. May panic if
/// unexpected behavior occurs.
pub(crate) fn mode_vec_usize(list: &Vec<usize>) -> Result<usize, Box<dyn Error>> {
    let mut counts: HashMap<usize, usize> = HashMap::new();

    for &num in list {
        *counts.entry(num).or_insert(1) += 1;
    }

    if counts.is_empty() {
        return Err(Box::from("Failed to calculate mode: input list is empty"));
    }

    Ok(*counts
        .iter()
        .max_by_key(|&(_, count)| count)
        .ok_or_else(|| Box::<dyn Error>::from("Failed to calculate mode: no mode found"))?
        .0)
}

#[allow(unused_imports)]
mod test {
    use super::*;
    use crate::fasta_manager::{open_fasta, Fasta};

    #[test]
    fn good_mode() {
        let the_list: Vec<usize> = Vec::from([2, 7, 9, 2, 7, 7, 3]);
        let mode = mode_vec_usize(&the_list);
        assert_eq!(mode.unwrap(), 7);
    }

    #[test]
    fn no_mode() {
        let the_list: Vec<usize> = Vec::new();
        let mode = mode_vec_usize(&the_list);
        assert_eq!(
            mode.unwrap_err().to_string(),
            "Failed to calculate mode: input list is empty"
        );
    }
}
