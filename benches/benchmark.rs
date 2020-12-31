use criterion::{black_box, criterion_group, criterion_main, Bencher, Criterion, Throughput};
use jigsaw::{
    AudioFrequency, F32SinSaw, InvMulSaw, IterativeSinSaw, Oscillator, ReferenceSaw, SamplingRateHz,
};

// Generate 1ms of oscillator signal, so that Criterion's per-second throughput
// measurement represents concurrent real-time signal generation abilities.
// (1k measurement / second = 1 real-time signal)
fn oscillator_benchmark<Osc: Oscillator>(
    sampling_rate: SamplingRateHz,
    freq: AudioFrequency,
) -> impl FnMut(&mut Bencher) {
    let mut oscillator = Osc::new(sampling_rate, freq, 0.0);
    let sample_count = sampling_rate / 1000 + (sampling_rate % 1000 != 0) as SamplingRateHz;
    move |b| {
        b.iter(|| {
            for _ in 0..sample_count {
                black_box(oscillator.next());
            }
        })
    }
}

// Benchmark sawtooth generators
pub fn saw_benchmark(criterion: &mut Criterion) {
    for sampling_rate in [44_100, 96_000, 192_000].iter().copied() {
        for saw_freq in [20.0, 200.0, 2_000.0, 20_000.0].iter().copied() {
            let group_title = format!("{} Hz saw, sampled at {} Hz", saw_freq, sampling_rate);
            let mut group = criterion.benchmark_group(group_title);
            group.throughput(Throughput::Elements(1));
            group.bench_function(
                "Reference",
                oscillator_benchmark::<ReferenceSaw>(sampling_rate, saw_freq),
            );
            group.bench_function(
                "f32 sinus",
                oscillator_benchmark::<F32SinSaw>(sampling_rate, saw_freq),
            );
            group.bench_function(
                "Iterative sinus",
                oscillator_benchmark::<IterativeSinSaw>(sampling_rate, saw_freq),
            );
            group.bench_function(
                "Multiply by inverse",
                oscillator_benchmark::<InvMulSaw>(sampling_rate, saw_freq),
            );
        }
    }
}

criterion_group!(benches, saw_benchmark);
criterion_main!(benches);
