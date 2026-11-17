//! Core logic for giraffe: indexing a GFF3 file in memory and querying it
//! with genome coordinates.

pub mod index;

pub use index::Index;

use std::error::Error;
use std::io::Read;

use csv::ReaderBuilder;
use multimap::MultiMap;
use serde::Deserialize;

pub type Result<T> = std::result::Result<T, Box<dyn Error>>;

/// A single query coordinate, read from a headerless `chrom<TAB>coord` table.
#[derive(Debug, Deserialize, Eq, PartialEq)]
pub struct Record {
    pub chrom: String,
    pub coord: u64,
}

/// Parse a headerless, tab-separated `chrom<TAB>coord` table.
pub fn read_coordinates<R: Read>(input: R) -> Result<MultiMap<String, u64>> {
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(false)
        .from_reader(input);

    let mut variants = MultiMap::new();

    for record in reader.deserialize() {
        let record: Record = record?;
        variants.insert(record.chrom, record.coord);
    }

    Ok(variants)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn read_coordinates_groups_by_chrom() {
        let v = read_coordinates("chr1\t10\nchr2\t20\nchr1\t30\n".as_bytes()).unwrap();
        assert_eq!(v.get_vec("chr1"), Some(&vec![10, 30]));
        assert_eq!(v.get_vec("chr2"), Some(&vec![20]));
    }

    #[test]
    fn read_coordinates_rejects_bad_input() {
        assert!(read_coordinates("chr1\tnotanumber\n".as_bytes()).is_err());
    }
}
