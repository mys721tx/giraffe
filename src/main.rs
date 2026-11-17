use std::fs::File;
use std::io::{self, BufReader, BufWriter, Read, Write};
use std::path::PathBuf;

use clap::{value_parser, Arg, Command};
use rusqlite::Connection;

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
                        .value_parser(value_parser!(PathBuf))
                        .action(clap::ArgAction::Set),
                )
                .arg(
                    Arg::new("input")
                        .short('i')
                        .long("input")
                        .value_name("IN")
                        .help("Path to the input tsv table. [default: stdin]")
                        .value_parser(value_parser!(PathBuf))
                        .action(clap::ArgAction::Set),
                )
                .arg(
                    Arg::new("output")
                        .short('o')
                        .long("output")
                        .value_name("OUT")
                        .help("Path to the output tsv table. [default: stdout]")
                        .value_parser(value_parser!(PathBuf))
                        .action(clap::ArgAction::Set),
                ),
        )
        .get_matches();

    match matches.subcommand() {
        Some(("build", matches)) => {
            let fin: Box<dyn Read> = match matches.get_one::<PathBuf>("input") {
                Some(f) => Box::new(BufReader::new(File::open(f).unwrap())),
                None => Box::new(io::stdin()),
            };

            let conn = Connection::open(matches.get_one::<PathBuf>("output").unwrap()).unwrap();

            giraffe::build(&conn, fin).unwrap();
        }
        Some(("query", matches)) => {
            let fin: Box<dyn Read> = match matches.get_one::<PathBuf>("input") {
                Some(f) => Box::new(BufReader::new(File::open(f).unwrap())),
                None => Box::new(io::stdin()),
            };

            let fout: Box<dyn Write> = match matches.get_one::<PathBuf>("output") {
                Some(f) => Box::new(BufWriter::new(File::create(f).unwrap())),
                None => Box::new(io::stdout()),
            };

            let conn = Connection::open(matches.get_one::<PathBuf>("database").unwrap()).unwrap();

            giraffe::query(&conn, fin, fout).unwrap();
        }
        _ => {}
    }
}
