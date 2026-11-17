use std::fs::File;
use std::io::{self, BufReader, BufWriter, Read, Write};
use std::path::PathBuf;

use clap::{value_parser, Arg, Command};
use giraffe::Index;

fn main() {
    let matches = Command::new("giraffe")
        .version("0.1.0")
        .author("Yishen Miao")
        .about("A GFF3 utility in Rust")
        .subcommand_required(true)
        .arg_required_else_help(true)
        .subcommand(
            Command::new("query")
                .about("Find the GFF3 features that overlap genome coordinates.")
                .arg(
                    Arg::new("gff")
                        .short('g')
                        .long("gff")
                        .value_name("GFF3")
                        .help("Path to the GFF3 annotation to search.")
                        .required(true)
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
                        .help("Path to the output GFF3 file. [default: stdout]")
                        .value_parser(value_parser!(PathBuf))
                        .action(clap::ArgAction::Set),
                ),
        )
        .get_matches();

    if let Some(("query", matches)) = matches.subcommand() {
        let fin: Box<dyn Read> = match matches.get_one::<PathBuf>("input") {
            Some(f) => Box::new(BufReader::new(File::open(f).unwrap())),
            None => Box::new(io::stdin()),
        };

        let fout: Box<dyn Write> = match matches.get_one::<PathBuf>("output") {
            Some(f) => Box::new(BufWriter::new(File::create(f).unwrap())),
            None => Box::new(BufWriter::new(io::stdout())),
        };

        let gff = matches.get_one::<PathBuf>("gff").unwrap();
        let index = Index::from_gff(BufReader::new(File::open(gff).unwrap())).unwrap();

        index.query(fin, fout).unwrap();
    }
}
