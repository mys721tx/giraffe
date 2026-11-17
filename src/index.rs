//! In-memory GFF3 index: records held in a `Vec`, with one interval tree per
//! sequence name pointing into it.

use std::collections::{BTreeMap, HashMap};
use std::io::{Read, Write};

use bio::data_structures::interval_tree::ArrayBackedIntervalTree;
use bio::io::gff;
use multimap::MultiMap;

use crate::Result;

/// GFF3 records indexed by sequence name and position.
///
/// GFF coordinates are 1-based and inclusive; the trees store them as the
/// half-open range `start..end + 1`.
#[derive(Default)]
pub struct Index {
    records: Vec<gff::Record>,
    trees: HashMap<String, ArrayBackedIntervalTree<u64, usize>>,
}

impl Index {
    /// Read every GFF3 record from `input` and index it.
    pub fn from_gff<R: Read>(input: R) -> Result<Self> {
        let mut index = Index::default();

        for r in gff::Reader::new(input, gff::GffType::GFF3).records() {
            let r = r?;
            let (start, end) = (*r.start(), *r.end());

            if end < start {
                return Err(format!(
                    "record {} on {} ends ({end}) before it starts ({start})",
                    index.records.len() + 1,
                    r.seqname()
                )
                .into());
            }

            index
                .trees
                .entry(r.seqname().to_owned())
                .or_default()
                .insert(start..end + 1, index.records.len());
            index.records.push(r);
        }

        index.trees.values_mut().for_each(|t| t.index());
        Ok(index)
    }

    /// Number of records in the index.
    pub fn len(&self) -> usize {
        self.records.len()
    }

    /// Whether the index holds no records.
    pub fn is_empty(&self) -> bool {
        self.records.is_empty()
    }

    /// The record with id `id` (its position in the input, from 0).
    pub fn record(&self, id: usize) -> Option<&gff::Record> {
        self.records.get(id)
    }

    /// Find every record overlapping each coordinate.
    ///
    /// Returns a map from record id to the coordinates that hit it, ordered
    /// by record id (that is, by input order).
    pub fn find_overlaps(&self, variants: &MultiMap<String, u64>) -> BTreeMap<usize, Vec<u64>> {
        let mut regions: BTreeMap<usize, Vec<u64>> = BTreeMap::new();
        // Reused across lookups; the entry type is not exported by `bio`.
        let mut buf = Vec::new();

        for (chrom, coords) in variants.iter_all() {
            let Some(tree) = self.trees.get(chrom) else {
                continue;
            };

            for &coord in coords {
                buf.clear();
                tree.find_into(coord..coord + 1, &mut buf);

                for e in &buf {
                    regions.entry(*e.data()).or_default().push(coord);
                }
            }
        }

        regions
    }

    /// Read coordinates from `input` and write every overlapping record to
    /// `output` as GFF3, in input order, with a `from` attribute listing the
    /// coordinates that hit it.
    ///
    /// Returns the number of records written.
    pub fn query<R: Read, W: Write>(&self, input: R, output: W) -> Result<usize> {
        let variants = crate::read_coordinates(input)?;
        let regions = self.find_overlaps(&variants);

        let mut writer = gff::Writer::new(output, gff::GffType::GFF3);

        for (&id, hits) in &regions {
            let mut r = self.records[id].clone();
            r.attributes_mut().insert(
                "from".to_string(),
                hits.iter()
                    .map(u64::to_string)
                    .collect::<Vec<_>>()
                    .join(" "),
            );
            writer.write(&r)?;
        }

        Ok(regions.len())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const GFF: &str = "\
##gff-version 3
chr1\tsrc\tgene\t100\t200\t.\t+\t.\tID=gene1;Name=Alpha
chr1\tsrc\texon\t150\t180\t5\t+\t0\tID=exon1;Parent=gene1
chr1\tsrc\tgene\t300\t400\t.\t-\t.\tID=gene2
chr2\tsrc\tgene\t100\t200\t.\t.\t.\tID=gene3
";

    fn index() -> Index {
        let index = Index::from_gff(GFF.as_bytes()).unwrap();
        assert_eq!(index.len(), 4);
        index
    }

    fn run_query(index: &Index, input: &str) -> Vec<gff::Record> {
        let mut out = Vec::new();
        index.query(input.as_bytes(), &mut out).unwrap();
        gff::Reader::new(out.as_slice(), gff::GffType::GFF3)
            .records()
            .map(|r| r.unwrap())
            .collect()
    }

    fn ids(records: &[gff::Record]) -> Vec<String> {
        records
            .iter()
            .map(|r| r.attributes().get("ID").unwrap().clone())
            .collect()
    }

    #[test]
    fn empty_input() {
        assert!(
            Index::from_gff("##gff-version 3\n".as_bytes())
                .unwrap()
                .is_empty()
        );
    }

    #[test]
    fn rejects_inverted_record() {
        let gff = "##gff-version 3\nchr1\tsrc\tgene\t200\t100\t.\t+\t.\tID=bad\n";
        assert!(Index::from_gff(gff.as_bytes()).is_err());
    }

    #[test]
    fn single_base_feature() {
        let gff = "##gff-version 3\nchr1\tsrc\tSNP\t50\t50\t.\t+\t.\tID=snp\n";
        let index = Index::from_gff(gff.as_bytes()).unwrap();
        assert_eq!(ids(&run_query(&index, "chr1\t50\n")), ["snp"]);
        assert!(run_query(&index, "chr1\t49\nchr1\t51\n").is_empty());
    }

    #[test]
    fn query_inside_nested_features() {
        let index = index();
        assert_eq!(ids(&run_query(&index, "chr1\t160\n")), ["gene1", "exon1"]);
    }

    #[test]
    fn query_bounds_are_inclusive() {
        let index = index();
        assert_eq!(ids(&run_query(&index, "chr1\t100\n")), ["gene1"]);
        assert_eq!(ids(&run_query(&index, "chr1\t200\n")), ["gene1"]);
        assert!(run_query(&index, "chr1\t99\n").is_empty());
        assert!(run_query(&index, "chr1\t201\n").is_empty());
    }

    #[test]
    fn query_respects_chromosome() {
        let index = index();
        assert_eq!(ids(&run_query(&index, "chr2\t150\n")), ["gene3"]);
        assert!(run_query(&index, "chr3\t150\n").is_empty());
    }

    #[test]
    fn query_no_hits_writes_nothing() {
        let index = index();
        let mut out = Vec::new();
        assert_eq!(index.query("chr1\t250\n".as_bytes(), &mut out).unwrap(), 0);
        assert!(out.is_empty());
    }

    #[test]
    fn query_records_from_attribute() {
        let index = index();
        let records = run_query(&index, "chr1\t110\nchr1\t160\nchr1\t350\n");

        let from: Vec<(String, String)> = records
            .iter()
            .map(|r| {
                (
                    r.attributes().get("ID").unwrap().clone(),
                    r.attributes().get("from").unwrap().clone(),
                )
            })
            .collect();

        assert_eq!(
            from,
            [
                ("gene1".to_string(), "110 160".to_string()),
                ("exon1".to_string(), "160".to_string()),
                ("gene2".to_string(), "350".to_string()),
            ]
        );
    }

    #[test]
    fn query_round_trips_fields() {
        let index = index();
        let records = run_query(&index, "chr1\t160\n");
        let exon = &records[1];

        assert_eq!(exon.seqname(), "chr1");
        assert_eq!(exon.source(), "src");
        assert_eq!(exon.feature_type(), "exon");
        assert_eq!((*exon.start(), *exon.end()), (150, 180));
        assert_eq!(exon.score(), Some(5));
        assert_eq!(format!("{:?}", exon.strand()), "Some(Forward)");
        assert_eq!(exon.phase(), &gff::Phase::try_from(0).unwrap());
        assert_eq!(exon.attributes().get("Parent").unwrap(), "gene1");

        let gene = &records[0];
        assert_eq!(gene.score(), None);
        assert_eq!(gene.phase(), &gff::Phase::try_from(None).unwrap());
        assert_eq!(gene.attributes().get("Name").unwrap(), "Alpha");
    }
}
