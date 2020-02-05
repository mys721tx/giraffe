use std::io;
use rusqlite::{named_params, Connection};
use bio::io::gff;

fn main() {

    let conn = Connection::open("anno.db").unwrap();

    conn.execute_batch(
        "CREATE TABLE IF NOT EXISTS anno (
            id INTEGER PRIMARY KEY,
            seqname TEXT,
            source TEXT,
            feature_type TEXT,
            start INTEGER,
            end INTEGER,
            score TEXT,
            strand TEXT,
            frame TEXT
        );
        CREATE TABLE IF NOT EXISTS attr (
            id INTEGER,
            type TEXT,
            value TEXT,
            PRIMARY KEY (id),
            FOREIGN KEY (id)
                REFERENCES anno (id)
                    ON DELETE CASCADE
                    ON UPDATE NO ACTION
        )").unwrap();

    let mut reader = gff::Reader::new(io::stdin(), gff::GffType::GFF3);

    let mut stmt = conn.prepare("INSERT INTO anno
        (seqname, source, feature_type, start, end, score, strand, frame)
        VALUES (:seqname, :source, :feature_type, :start, :end, :score, :strand, :frame)").unwrap();

    conn.execute_batch("BEGIN TRANSACTION;").unwrap();

    for r in reader.records() {
        let r = r.unwrap();

        let start = *r.start() as i64;
        let end = *r.end() as i64;
        let score = r.score().map(|x| x as i64);
        let strand = r.strand().map(|x| x.to_string());

        stmt.execute_named(
            named_params!{
                ":seqname": r.seqname(),
                ":source": r.source(),
                ":feature_type": r.feature_type(),
                ":start": start,
                ":end": end,
                ":score": score,
                ":strand": strand,
                ":frame": r.frame()
            }
        ).unwrap();
    }

    conn.execute_batch("END TRANSACTION;").unwrap();
}
