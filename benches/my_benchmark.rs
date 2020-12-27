use criterion::{criterion_group, criterion_main, Criterion, Throughput};
use jigsaw::ReferenceSaw;

pub fn saw_benchmark(criterion: &mut Criterion) {
    for sample_rate in [44_100, 96_000, 192_000].iter().copied() {
        for saw_freq in [20.0, 200.0, 2_000.0, 20_000.0].iter().copied() {
            let group_title = format!("{} Hz saw, sampled at {} Hz", saw_freq, sample_rate);
            let mut group = criterion.benchmark_group(group_title);
            group.throughput(Throughput::Elements(1));

            let mut reference_saw = ReferenceSaw::new(sample_rate, saw_freq, 0.0);
            group.bench_function("Reference", |b| b.iter(|| reference_saw.next()));
        }
    }
}

criterion_group!(benches, saw_benchmark);
criterion_main!(benches);
