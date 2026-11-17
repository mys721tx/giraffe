use std::fmt::Write as _;
use std::hint::black_box;

use criterion::{criterion_group, criterion_main, BatchSize, BenchmarkId, Criterion, Throughput};
use rusqlite::Connection;

const CHROMS: usize = 4;

/// Generate `n` genes, each with two exons, spread over `CHROMS` chromosomes.
fn synthetic_gff(n: usize) -> String {
    let mut s = String::from("##gff-version 3\n");
    for i in 0..n {
        let chrom = i % CHROMS;
        let start = (i / CHROMS) * 1_000 + 1;
        let end = start + 800;
        let strand = if i % 2 == 0 { '+' } else { '-' };
        writeln!(
            s,
            "chr{chrom}\tbench\tgene\t{start}\t{end}\t.\t{strand}\t.\tID=gene{i};Name=G{i}"
        )
        .unwrap();
        for (e, (es, ee)) in [(start, start + 199), (start + 500, end)]
            .iter()
            .enumerate()
        {
            writeln!(
                s,
                "chr{chrom}\tbench\texon\t{es}\t{ee}\t.\t{strand}\t0\tID=exon{i}_{e};Parent=gene{i}"
            )
            .unwrap();
        }
    }
    s
}

/// Generate `n` query coordinates: alternating hits (inside an exon) and misses.
fn synthetic_coords(n: usize, genes: usize) -> String {
    let mut s = String::new();
    for i in 0..n {
        let g = (i * 7919) % genes;
        let chrom = g % CHROMS;
        let base = (g / CHROMS) * 1_000 + 1;
        let coord = if i % 2 == 0 { base + 100 } else { base + 900 };
        writeln!(s, "chr{chrom}\t{coord}").unwrap();
    }
    s
}

fn build_db(gff: &str) -> Connection {
    let conn = Connection::open_in_memory().unwrap();
    giraffe::build(&conn, gff.as_bytes()).unwrap();
    conn
}

fn bench_build(c: &mut Criterion) {
    let mut group = c.benchmark_group("build");
    group.sample_size(20);

    for genes in [1_000, 10_000] {
        let gff = synthetic_gff(genes);
        group.throughput(Throughput::Elements((genes * 3) as u64));
        group.bench_with_input(BenchmarkId::from_parameter(genes * 3), &gff, |b, gff| {
            b.iter_batched(
                || Connection::open_in_memory().unwrap(),
                |conn| giraffe::build(&conn, black_box(gff.as_bytes())).unwrap(),
                BatchSize::LargeInput,
            )
        });
    }

    group.finish();
}

fn bench_query(c: &mut Criterion) {
    const GENES: usize = 10_000;
    let conn = build_db(&synthetic_gff(GENES));

    let mut group = c.benchmark_group("query");

    for n in [100, 1_000, 10_000] {
        let coords = synthetic_coords(n, GENES);
        group.throughput(Throughput::Elements(n as u64));
        group.bench_with_input(BenchmarkId::from_parameter(n), &coords, |b, coords| {
            let mut out = Vec::with_capacity(1 << 20);
            b.iter(|| {
                out.clear();
                giraffe::query(&conn, black_box(coords.as_bytes()), &mut out).unwrap()
            })
        });
    }

    group.finish();
}

fn bench_find_overlaps(c: &mut Criterion) {
    const GENES: usize = 10_000;
    let conn = build_db(&synthetic_gff(GENES));
    let variants = giraffe::read_coordinates(synthetic_coords(1_000, GENES).as_bytes()).unwrap();

    c.bench_function("find_overlaps/1000", |b| {
        b.iter(|| giraffe::find_overlaps(&conn, black_box(&variants)).unwrap())
    });
}

criterion_group!(benches, bench_build, bench_query, bench_find_overlaps);
criterion_main!(benches);
