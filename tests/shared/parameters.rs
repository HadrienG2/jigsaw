//! All integration tests work in the same parameter space

use core::ops::Range;
use genawaiter::yield_;
use jigsaw::{AudioFrequency, AudioPhase, AudioPhaseMod, SamplingRateHz};
use rand::Rng;

// Region of interest within the oscillator configuration space
pub const SAMPLING_RATES: &'static [SamplingRateHz] = &[44_100, 48_000, 96_000, 192_000];
pub const OSCILLATOR_FREQ_RANGE: Range<AudioFrequency> = 20.0..20000.0;
pub const PHASE_RANGE: Range<AudioPhase> = 0.0..AudioPhaseMod::consts::TAU;

/// Range of the base-2 logarithm of the relative sample rate
///
/// The relative sample rate is proportional to the number of harmonics in the
/// band-limited sawtooth wave, and the error on a band-limited approximation of
/// an infinite Fourier series converges logarithmically as the number of
/// harmonics in said Fourier series increases. So to observe the convergence of
/// the Fourier series, it is best to sample the relative frequency range
/// logarithmically.
///
/// Among all possible logarithm bases, base-2 logarithms should be computable
/// with optimal precision for IEEE-754 binary floats.
///
pub fn log2_relative_rate_range() -> Range<AudioFrequency> {
    let start =
        ((*SAMPLING_RATES.first().unwrap() as AudioFrequency) / OSCILLATOR_FREQ_RANGE.end).log2();
    let end =
        ((*SAMPLING_RATES.last().unwrap() as AudioFrequency) / OSCILLATOR_FREQ_RANGE.start).log2();
    start..end
}

// Number of phase and relative frequency buckets
//
// This gives the spatial resolution of the final map of the sawtooth wave's
// properties. More buckets take more time to compute, but allow studying finer
// details of the parameter space.
//
pub const NUM_PHASE_BUCKETS: usize = 1000;
pub const NUM_RELATIVE_FREQ_BUCKETS: usize = 2000;

/// Size of a certain parameter's buckets
fn bucket_size(range: Range<f32>, num_buckets: usize) -> f32 {
    (range.end - range.start) / (num_buckets as f32)
}

/// Start of a certain bucket within a parameter space
pub fn bucket_start(range: Range<f32>, num_buckets: usize, idx: usize) -> f32 {
    range.start + (idx as f32) * bucket_size(range, num_buckets)
}

/// Number of data points to be measured in every bucket
///
/// More samples samples per bucket will produce a less noisy map of the
/// sawtooth wave's properties, in a manner similar to how super-sampling
/// antialiasing reduces aliasing artifacts in 3D graphics. As with
/// super-sampling antialiasing, the price to pay is a quadratic increase of
/// a test's running time.
///
pub const SAMPLES_PER_BUCKET: usize = 2;

/// Given a certain parameter range, which is divided into a certain number of
/// evenly sized buckets, produce a set of SAMPLES_PER_BUCKET random parameter
/// values falling inside of each bucket.
///
/// The associated bucket index is provided for convenience.
///
pub fn irregular_samples(
    range: Range<f32>,
    num_buckets: usize,
) -> impl Iterator<Item = (usize, [f32; SAMPLES_PER_BUCKET])> {
    genawaiter::rc::gen!({
        let mut rng = rand::thread_rng();
        let bucket_size = bucket_size(range.clone(), num_buckets);
        let mut samples = [0.0; SAMPLES_PER_BUCKET];
        for bucket in 0..num_buckets {
            let bucket_start = bucket_start(range.clone(), num_buckets, bucket);
            let bucket_end = bucket_start + bucket_size;
            for sample in &mut samples[..] {
                *sample = rng.gen_range(bucket_start..bucket_end);
            }
            yield_!((bucket, samples));
        }
    })
    .into_iter()
}
