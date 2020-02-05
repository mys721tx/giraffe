use std::io;
use rusqlite::{params, Connection};
use bio::io::gff;

fn main() {

    let conn = Connection::open("anno.db").unwrap();

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
        );
        ").unwrap();

    let mut reader = gff::Reader::new(io::stdin(), gff::GffType::GFF3);

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
}
