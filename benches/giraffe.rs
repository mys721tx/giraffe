use std::fmt::Write as _;
use std::hint::black_box;

use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion, Throughput};
use giraffe::Index;

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

fn bench_build(c: &mut Criterion) {
    let mut group = c.benchmark_group("build");
    group.sample_size(20);

    for genes in [1_000, 10_000] {
        let gff = synthetic_gff(genes);
        group.throughput(Throughput::Elements((genes * 3) as u64));
        group.bench_with_input(BenchmarkId::from_parameter(genes * 3), &gff, |b, gff| {
            b.iter(|| Index::from_gff(black_box(gff.as_bytes())).unwrap())
        });
    }

    group.finish();
}

fn bench_query(c: &mut Criterion) {
    const GENES: usize = 10_000;
    let index = Index::from_gff(synthetic_gff(GENES).as_bytes()).unwrap();

    let mut group = c.benchmark_group("query");

    for n in [100, 1_000, 10_000] {
        let coords = synthetic_coords(n, GENES);
        group.throughput(Throughput::Elements(n as u64));
        group.bench_with_input(BenchmarkId::from_parameter(n), &coords, |b, coords| {
            let mut out = Vec::with_capacity(1 << 20);
            b.iter(|| {
                out.clear();
                index.query(black_box(coords.as_bytes()), &mut out).unwrap()
            })
        });
    }

    group.finish();
}

fn bench_find_overlaps(c: &mut Criterion) {
    const GENES: usize = 10_000;
    let index = Index::from_gff(synthetic_gff(GENES).as_bytes()).unwrap();
    let variants = giraffe::read_coordinates(synthetic_coords(1_000, GENES).as_bytes()).unwrap();

    c.bench_function("find_overlaps/1000", |b| {
        b.iter(|| index.find_overlaps(black_box(&variants)))
    });
}

criterion_group!(benches, bench_build, bench_query, bench_find_overlaps);
criterion_main!(benches);
