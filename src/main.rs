use std::fs::File;
use std::io;
use std::io::{BufReader, BufWriter, Read, Write};
use std::mem::replace;
use std::path::PathBuf;

use bio::io::gff;
use clap::{value_parser, Arg, Command};
use csv::ReaderBuilder;
use multimap::MultiMap;
use rusqlite::types::Null;
use rusqlite::{params, Connection};
use serde::Deserialize;

#[derive(Debug, Deserialize, Eq, PartialEq)]
struct Record {
    chrom: String,
    coord: u64,
}

fn main() {
    let matches = Command::new("giraffe")
        .version("0.1.0")
        .author("Yishen Miao")
        .about("A GFF3 utility in Rust")
        .subcommand(
            Command::new("build")
                .about("Build a SQLite database from a GFF3 file.")
                .arg(
                    Arg::new("input")
                        .short('i')
                        .long("input")
                        .value_name("GFF3")
                        .value_parser(value_parser!(PathBuf))
                        .help("Sets the name of the input file. Defaults to stdin.")
                        .action(clap::ArgAction::Set),
                )
                .arg(
                    Arg::new("output")
                        .short('o')
                        .long("output")
                        .value_name("DB")
                        .help("Sets the name of the output file. Defaults to \"anno.db\".")
                        .default_value("anno.db")
                        .value_parser(value_parser!(PathBuf))
                        .action(clap::ArgAction::Set),
                ),
        )
        .subcommand(
            Command::new("query")
                .about("Query genome coordinates in a SQLite database.")
                .arg(
                    Arg::new("database")
                        .short('d')
                        .long("db")
                        .value_name("DB")
                        .help("Path to the SQLite database.")
                        .default_value("anno.db")
                        .action(clap::ArgAction::Set),
                )
                .arg(
                    Arg::new("input")
                        .short('i')
                        .long("input")
                        .value_name("IN")
                        .help("Path to the input tsv table. [default: stdin]")
                        .action(clap::ArgAction::Set),
                )
                .arg(
                    Arg::new("output")
                        .short('o')
                        .long("output")
                        .value_name("OUT")
                        .help("Path to the output tsv table. [default: stdout]")
                        .action(clap::ArgAction::Set),
                ),
        )
        .get_matches();

    match matches.subcommand_name() {
        Some("build") => {
            let matches = matches.subcommand_matches("build").unwrap();

            let fin: Box<dyn Read> = match matches.get_one::<PathBuf>("input") {
                Some(f) => Box::new(BufReader::new(File::open(f).unwrap())),
                None => Box::new(io::stdin()),
            };

            let conn = Connection::open(matches.get_one::<PathBuf>("output").unwrap()).unwrap();

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
                    +frame TEXT
                );
                CREATE TABLE IF NOT EXISTS attr (
                    anno_id INTEGER,
                    type TEXT,
                    value TEXT,
                    FOREIGN KEY (anno_id)
                        REFERENCES anno (id)
                            ON DELETE CASCADE
                            ON UPDATE NO ACTION
                );
                CREATE INDEX anno_id ON attr (anno_id);",
            )
            .unwrap();

            let mut reader = gff::Reader::new(fin, gff::GffType::GFF3);

            let mut stmt_anno = conn
                .prepare(
                    "INSERT INTO anno (
                        id,
                        start,
                        end,
                        seqname,
                        source,
                        feature_type,
                        score,
                        strand,
                        frame
                    ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)",
                )
                .unwrap();

            let mut stmt_attr = conn
                .prepare("INSERT INTO attr (anno_id, type, value) VALUES (?, ?, ?)")
                .unwrap();

            conn.execute_batch("BEGIN TRANSACTION;").unwrap();

            for r in reader.records() {
                let r = r.unwrap();

                let start = *r.start() as i64;
                let end = *r.end() as i64;
                let score = r.score().map(|x| x as i64);
                let strand = r.strand().map(|x| x.to_string());

                let row = stmt_anno
                    .insert(params![
                        Null,
                        start,
                        end,
                        r.seqname(),
                        r.source(),
                        r.feature_type(),
                        score,
                        strand,
                        r.frame()
                    ])
                    .unwrap();

                for (key, value) in r.attributes().iter() {
                    stmt_attr.insert(params![row, key, value]).unwrap();
                }
            }

            conn.execute_batch("END TRANSACTION;").unwrap();
        }
        Some("query") => {
            let matches = matches.subcommand_matches("query").unwrap();

            let fin: Box<dyn Read> = match matches.get_one::<PathBuf>("input") {
                Some(f) => Box::new(BufReader::new(File::open(f).unwrap())),
                None => Box::new(io::stdin()),
            };

            let fout: Box<dyn Write> = match matches.get_one::<PathBuf>("output") {
                Some(f) => Box::new(BufWriter::new(File::create(f).unwrap())),
                None => Box::new(io::stdout()),
            };

            let conn = Connection::open(matches.get_one::<PathBuf>("database").unwrap()).unwrap();

            let mut reader = ReaderBuilder::new()
                .delimiter(b'\t')
                .has_headers(false)
                .from_reader(fin);

            let mut variants = MultiMap::new();

            for record in reader.deserialize() {
                let record: Record = record.unwrap();
                variants.insert(record.chrom, record.coord);
            }

            let mut stmt_get_intervals = conn
                .prepare(
                    "SELECT id FROM anno
                    WHERE seqname = (?)
                        AND start <= (?)
                        AND end >= (?)
                ",
                )
                .unwrap();

            let mut regions = MultiMap::new();

            for (chrom, coords) in variants.iter_all() {
                for coord in coords.iter() {
                    let coord = *coord as i64;

                    let mut rows = stmt_get_intervals
                        .query(params![chrom, coord, coord])
                        .unwrap();

                    while let Some(row) = rows.next().unwrap() {
                        let anno_id: i64 = row.get_unwrap(0);
                        regions.insert(anno_id, coord);
                    }
                }
            }

            let mut stmt_get_anno = conn.prepare(
                "SELECT seqname, source, feature_type, score, start, end, strand, frame FROM anno
                    WHERE id = (?)"
            ).unwrap();

            let mut stmt_get_attr = conn
                .prepare("SELECT type, value FROM attr WHERE anno_id = (?)")
                .unwrap();

            let mut writer = gff::Writer::new(fout, gff::GffType::GFF3);

            for id in regions.keys() {
                let mut rows = stmt_get_anno.query(params![id]).unwrap();

                while let Some(row) = rows.next().unwrap() {
                    let mut r = gff::Record::new();

                    let _ = replace(r.seqname_mut(), row.get_unwrap::<usize, String>(0));

                    let _ = replace(r.source_mut(), row.get_unwrap::<usize, String>(1));

                    let _ = replace(r.feature_type_mut(), row.get_unwrap::<usize, String>(2));

                    let _ = replace(
                        r.score_mut(),
                        row.get::<usize, i64>(3)
                            .ok()
                            .map_or_else(|| ".".to_string(), |x| x.to_string()),
                    );

                    let _ = replace(r.start_mut(), row.get_unwrap::<usize, i64>(4) as u64);

                    let _ = replace(r.end_mut(), row.get_unwrap::<usize, i64>(5) as u64);

                    let _ = replace(
                        r.strand_mut(),
                        row.get::<usize, String>(6)
                            .ok()
                            .unwrap_or_else(|| ".".to_string()),
                    );

                    let _ = replace(
                        r.frame_mut(),
                        row.get::<usize, String>(8)
                            .ok()
                            .unwrap_or_else(|| "".to_string()),
                    );

                    let mut attributes: MultiMap<String, String> = MultiMap::new();

                    let mut rows = stmt_get_attr.query(params![id]).unwrap();

                    while let Some(row) = rows.next().unwrap() {
                        let key = row.get_unwrap(0);
                        let value = row.get_unwrap(1);

                        attributes.insert(key, value);
                    }

                    attributes.insert(
                        "from".to_string(),
                        regions
                            .get_vec(id)
                            .unwrap()
                            .iter()
                            .map(|x| x.to_string())
                            .collect::<Vec<String>>()
                            .join(" "),
                    );

                    let _ = replace(r.attributes_mut(), attributes);

                    writer.write(&r).unwrap();
                }
            }
        }
        _ => {}
    }
}
