use std::io;
use std::fs::File;
use std::io::Read;

use clap::{Arg, App, SubCommand};
use rusqlite::{params, Connection};
use bio::io::gff;

fn main() {

    let matches = App::new("giraffe")
        .version("0.1.0")
        .author("Yishen Miao")
        .about("A GFF3 utility in Rust")
        .subcommand(SubCommand::with_name("build")
            .about("Build a SQLite database from a GFF3 file.")
            .arg(Arg::with_name("input")
                .short("i")
                .long("input")
                .value_name("GFF3")
                .help("Path to the GFF3 file. [default: stdin]")
                .takes_value(true)
            )
            .arg(Arg::with_name("output")
                .short("o")
                .long("output")
                .value_name("DB")
                .help("Path to the SQLite database.")
                .default_value("anno.db")
                .takes_value(true)
            )
        )
        .get_matches();

    match matches.subcommand_name() {
        Some("build")  => {

            let matches = matches.subcommand_matches("build").unwrap();

            let fin: Box<dyn Read> = match matches.value_of("input") {
                Some(f) => Box::new(File::open(f).unwrap()),
                None => Box::new(io::stdin()),
            };

            let conn = Connection::open(matches.value_of("output").unwrap()).unwrap();

            conn.execute_batch(
                "CREATE TABLE IF NOT EXISTS anno (
                    id INTEGER PRIMARY KEY,
                    seqname TEXT,
                    source TEXT,
                    feature_type TEXT,
                    score TEXT,
                    strand TEXT,
                    frame TEXT
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
                CREATE INDEX anno_id ON attr (anno_id);
                CREATE VIRTUAL TABLE domain USING rtree_i32(
                    id,
                    start,
                    end
                );").unwrap();

            let mut reader = gff::Reader::new(fin, gff::GffType::GFF3);

            let mut stmt_anno = conn.prepare("INSERT INTO anno
                (seqname, source, feature_type, score, strand, frame)
                VALUES (?, ?, ?, ?, ?, ?)").unwrap();

            let mut stmt_attr = conn.prepare("INSERT INTO attr
                (anno_id, type, value) VALUES (?, ?, ?)").unwrap();

            let mut stmt_rtree = conn.prepare("INSERT INTO domain
                (id, start, end)
                VALUES (?, ?, ?)").unwrap();

            conn.execute_batch("BEGIN TRANSACTION;").unwrap();

            for r in reader.records() {
                let r = r.unwrap();

                let start = *r.start() as i64;
                let end = *r.end() as i64;
                let score = r.score().map(|x| x as i64);
                let strand = r.strand().map(|x| x.to_string());

                let row = stmt_anno.insert(
                    params![
                        r.seqname(),
                        r.source(),
                        r.feature_type(),
                        score,
                        strand,
                        r.frame()
                    ]
                ).unwrap();

                stmt_rtree.insert(params![row, start, end]).unwrap();

                for (key, value) in r.attributes().iter() {
                    stmt_attr.insert(params![row, key, value]).unwrap();
                }
            }

            conn.execute_batch("END TRANSACTION;").unwrap();
        },
        Some("query") => {},
        _ => {},
    }
}
