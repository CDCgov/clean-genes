#![allow(clippy::should_panic_without_expect)]
use std::collections::HashMap;

/// Calculates the mathematical mode of a vector of usizes.
pub(crate) fn mode_vec_usize(list: &Vec<usize>) -> Option<usize> {
    let mut counts: HashMap<usize, usize> = HashMap::new();

    for &num in list {
        *counts.entry(num).or_default() += 1;
    }

    if counts.is_empty() {
        return None;
    }

    let mode = *counts.iter().max_by_key(|&(_, count)| count)?.0;
    Some(mode)
}

#[expect(unused_imports)]
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
    #[should_panic]
    fn no_mode() {
        let the_list: Vec<usize> = Vec::new();
        let mode = mode_vec_usize(&the_list);
        mode.expect("Failed to calculate mode: input list is empty");
    }
}
