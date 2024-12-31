use std::collections::HashMap;

/// Calculates the mathematical mode of a vector of usizes
pub(crate) fn get_mode_vec_usize(list: &Vec<usize>) -> Option<usize> {
    let mut counts: HashMap<usize, usize> = HashMap::new();

    for &num in list {
        if let Some(value) = counts.get_mut(&num) {
            *value += 1;
        } else {
            counts.insert(num, 1);
        }
    }

    if counts.is_empty() {
        None
    } else {
        counts
            .into_iter()
            .max_by_key(|&(_, count)| count)
            .map(|(number, _)| number)
    }
}

#[allow(unused_imports)]
mod test {
    use super::*;
    use crate::fasta_manager::{open_fasta, Fasta};

    #[test]
    fn good_mode() {
        let the_list: Vec<usize> = vec![2, 7, 9, 2, 7, 7, 3];
        let mode = get_mode_vec_usize(&the_list);
        assert_eq!(mode, Some(7));
    }

    #[test]
    fn no_mode() {
        let the_list: Vec<usize> = Vec::new();
        let mode = get_mode_vec_usize(&the_list);
        assert_eq!(mode, None);
    }
}
