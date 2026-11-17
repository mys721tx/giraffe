//! Core logic for giraffe: building an R*Tree-indexed SQLite database from a
//! GFF3 file and querying it with genome coordinates.

use std::error::Error;
use std::io::{Read, Write};

use bio::io::gff;
use csv::ReaderBuilder;
use multimap::MultiMap;
use rusqlite::types::Null;
use rusqlite::{params, Connection};
use serde::Deserialize;

pub type Result<T> = std::result::Result<T, Box<dyn Error>>;

/// A single query coordinate, read from a headerless `chrom<TAB>coord` table.
#[derive(Debug, Deserialize, Eq, PartialEq)]
pub struct Record {
    pub chrom: String,
    pub coord: u64,
}

/// Create the `anno` R*Tree table and the `attr` attribute table.
pub fn create_schema(conn: &Connection) -> Result<()> {
    conn.execute_batch(
        "CREATE VIRTUAL TABLE anno USING rtree_i32 (
            id,
            start,
            end,
            +seqname TEXT,
            +source TEXT,
            +feature_type TEXT,
            +score TEXT,
            +strand TEXT,
            +phase
        );
        CREATE TABLE IF NOT EXISTS attr (
            anno_id INTEGER,
            type TEXT,
            value TEXT
        );
        CREATE INDEX anno_id ON attr (anno_id);",
    )?;
    Ok(())
}

/// Create the schema in `conn` and load every GFF3 record from `input`.
///
/// Returns the number of records inserted.
pub fn build<R: Read>(conn: &Connection, input: R) -> Result<usize> {
    create_schema(conn)?;

    let tx = conn.unchecked_transaction()?;
    let mut count = 0;

    {
        let mut stmt_anno = tx.prepare(
            "INSERT INTO anno (
                id,
                start,
                end,
                seqname,
                source,
                feature_type,
                score,
                strand,
                phase
            ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)",
        )?;

        let mut stmt_attr =
            tx.prepare("INSERT INTO attr (anno_id, type, value) VALUES (?, ?, ?)")?;

        let mut reader = gff::Reader::new(input, gff::GffType::GFF3);

        for r in reader.records() {
            let r = r?;

            let start = *r.start() as i64;
            let end = *r.end() as i64;
            let score = r.score().map(|x| x as i64);
            let strand = r.strand().map(|x| x.to_string());
            let phase: Option<u8> = r.phase().clone().try_into().ok();

            let row = stmt_anno.insert(params![
                Null,
                start,
                end,
                r.seqname(),
                r.source(),
                r.feature_type(),
                score,
                strand,
                phase
            ])?;

            for (key, value) in r.attributes().iter() {
                stmt_attr.insert(params![row, key, value])?;
            }

            count += 1;
        }
    }

    tx.commit()?;
    Ok(count)
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

/// Find every annotation overlapping each coordinate.
///
/// Returns a map from annotation id to the coordinates that hit it.
pub fn find_overlaps(
    conn: &Connection,
    variants: &MultiMap<String, u64>,
) -> Result<MultiMap<i64, i64>> {
    let mut stmt = conn.prepare_cached(
        "SELECT id FROM anno
        WHERE seqname = (?)
            AND start <= (?)
            AND end >= (?)",
    )?;

    let mut regions = MultiMap::new();

    for (chrom, coords) in variants.iter_all() {
        for coord in coords.iter() {
            let coord = *coord as i64;

            let mut rows = stmt.query(params![chrom, coord, coord])?;

            while let Some(row) = rows.next()? {
                let anno_id: i64 = row.get(0)?;
                regions.insert(anno_id, coord);
            }
        }
    }

    Ok(regions)
}

/// Load the annotation `id` as a GFF record, adding a `from` attribute that
/// lists the query coordinates `hits` that overlapped it.
pub fn fetch_record(conn: &Connection, id: i64, hits: &[i64]) -> Result<Option<gff::Record>> {
    let mut stmt_get_anno = conn.prepare_cached(
        "SELECT seqname, source, feature_type, score, start, end, strand, phase FROM anno
        WHERE id = (?)",
    )?;

    let mut stmt_get_attr =
        conn.prepare_cached("SELECT type, value FROM attr WHERE anno_id = (?)")?;

    let mut rows = stmt_get_anno.query(params![id])?;

    let row = match rows.next()? {
        Some(row) => row,
        None => return Ok(None),
    };

    let mut r = gff::Record::new();

    *r.seqname_mut() = row.get(0)?;
    *r.source_mut() = row.get(1)?;
    *r.feature_type_mut() = row.get(2)?;
    *r.score_mut() = row
        .get::<_, i64>(3)
        .ok()
        .map_or_else(|| ".".to_string(), |x| x.to_string());
    *r.start_mut() = row.get::<_, i64>(4)? as u64;
    *r.end_mut() = row.get::<_, i64>(5)? as u64;
    *r.strand_mut() = row
        .get::<_, String>(6)
        .ok()
        .unwrap_or_else(|| ".".to_string());
    *r.phase_mut() = gff::Phase::try_from(row.get::<_, Option<u8>>(7)?)?;

    let mut attributes: MultiMap<String, String> = MultiMap::new();

    let mut attr_rows = stmt_get_attr.query(params![id])?;

    while let Some(row) = attr_rows.next()? {
        attributes.insert(row.get(0)?, row.get(1)?);
    }

    attributes.insert(
        "from".to_string(),
        hits.iter()
            .map(|x| x.to_string())
            .collect::<Vec<String>>()
            .join(" "),
    );

    *r.attributes_mut() = attributes;

    Ok(Some(r))
}

/// Read coordinates from `input`, look them up in `conn` and write every
/// overlapping annotation to `output` as GFF3, ordered by annotation id.
///
/// Returns the number of records written.
pub fn query<R: Read, W: Write>(conn: &Connection, input: R, output: W) -> Result<usize> {
    let variants = read_coordinates(input)?;
    let regions = find_overlaps(conn, &variants)?;

    let mut ids: Vec<i64> = regions.keys().copied().collect();
    ids.sort_unstable();

    let mut writer = gff::Writer::new(output, gff::GffType::GFF3);
    let mut count = 0;

    for id in ids {
        let hits = regions.get_vec(&id).map(Vec::as_slice).unwrap_or(&[]);

        if let Some(r) = fetch_record(conn, id, hits)? {
            writer.write(&r)?;
            count += 1;
        }
    }

    Ok(count)
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

    fn db() -> Connection {
        let conn = Connection::open_in_memory().unwrap();
        assert_eq!(build(&conn, GFF.as_bytes()).unwrap(), 4);
        conn
    }

    fn run_query(conn: &Connection, input: &str) -> Vec<gff::Record> {
        let mut out = Vec::new();
        query(conn, input.as_bytes(), &mut out).unwrap();
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
    fn build_stores_records() {
        let conn = db();

        let n: i64 = conn
            .query_row("SELECT count(*) FROM anno", [], |r| r.get(0))
            .unwrap();
        assert_eq!(n, 4);

        let (start, end, seqname, score, strand, phase): (
            i64,
            i64,
            String,
            Option<i64>,
            Option<String>,
            Option<u8>,
        ) = conn
            .query_row(
                "SELECT a.start, a.end, a.seqname, a.score, a.strand, a.phase
                FROM anno a JOIN attr t ON t.anno_id = a.id
                WHERE t.type = 'ID' AND t.value = 'exon1'",
                [],
                |r| {
                    Ok((
                        r.get(0)?,
                        r.get(1)?,
                        r.get(2)?,
                        r.get(3)?,
                        r.get(4)?,
                        r.get(5)?,
                    ))
                },
            )
            .unwrap();

        assert_eq!((start, end), (150, 180));
        assert_eq!(seqname, "chr1");
        assert_eq!(score, Some(5));
        assert_eq!(strand.as_deref(), Some("+"));
        assert_eq!(phase, Some(0));
    }

    #[test]
    fn build_stores_attributes() {
        let conn = db();

        let n: i64 = conn
            .query_row("SELECT count(*) FROM attr", [], |r| r.get(0))
            .unwrap();
        // gene1: ID, Name; exon1: ID, Parent; gene2: ID; gene3: ID
        assert_eq!(n, 6);
    }

    #[test]
    fn build_unknown_strand_is_null() {
        let conn = db();

        let strand: Option<String> = conn
            .query_row(
                "SELECT a.strand FROM anno a JOIN attr t ON t.anno_id = a.id
                WHERE t.value = 'gene3'",
                [],
                |r| r.get(0),
            )
            .unwrap();
        assert_eq!(strand, None);
    }

    #[test]
    fn build_twice_fails() {
        let conn = db();
        assert!(build(&conn, GFF.as_bytes()).is_err());
    }

    #[test]
    fn build_empty_input() {
        let conn = Connection::open_in_memory().unwrap();
        assert_eq!(build(&conn, "##gff-version 3\n".as_bytes()).unwrap(), 0);
    }

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

    #[test]
    fn query_inside_nested_features() {
        let conn = db();
        let records = run_query(&conn, "chr1\t160\n");
        assert_eq!(ids(&records), ["gene1", "exon1"]);
    }

    #[test]
    fn query_bounds_are_inclusive() {
        let conn = db();
        assert_eq!(ids(&run_query(&conn, "chr1\t100\n")), ["gene1"]);
        assert_eq!(ids(&run_query(&conn, "chr1\t200\n")), ["gene1"]);
        assert!(run_query(&conn, "chr1\t99\n").is_empty());
        assert!(run_query(&conn, "chr1\t201\n").is_empty());
    }

    #[test]
    fn query_respects_chromosome() {
        let conn = db();
        assert_eq!(ids(&run_query(&conn, "chr2\t150\n")), ["gene3"]);
        assert!(run_query(&conn, "chr3\t150\n").is_empty());
    }

    #[test]
    fn query_no_hits_writes_nothing() {
        let conn = db();
        let mut out = Vec::new();
        assert_eq!(query(&conn, "chr1\t250\n".as_bytes(), &mut out).unwrap(), 0);
        assert!(out.is_empty());
    }

    #[test]
    fn query_records_from_attribute() {
        let conn = db();
        let records = run_query(&conn, "chr1\t110\nchr1\t160\nchr1\t350\n");

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
        let conn = db();
        let records = run_query(&conn, "chr1\t160\n");
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

    #[test]
    fn fetch_record_missing_id() {
        let conn = db();
        assert!(fetch_record(&conn, 9999, &[]).unwrap().is_none());
    }
}
