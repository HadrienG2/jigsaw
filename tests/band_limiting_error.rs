//! This test studies the error of band-limited signals, both with respect to
//! each other and to band-unlimited signals.

use core::sync::atomic::{AtomicU8, Ordering};
use genawaiter::yield_;
use jigsaw::{
    unlimited_saw, AudioFrequency, AudioPhase, AudioPhaseMod, AudioSample, Oscillator,
    ReferenceSaw, SamplingRateHz,
};
use log::{debug, trace};
use plotters::prelude::*;
use rand::Rng;
use rayon::prelude::*;
use std::ops::Range;

// Region of interest within the oscillator configuration space
const SAMPLING_RATE_RANGE: Range<SamplingRateHz> = 44_100..192_000;
const OSCILLATOR_FREQ_RANGE: Range<AudioFrequency> = 20.0..20000.0;
const PHASE_RANGE: Range<AudioPhase> = 0.0..AudioPhaseMod::consts::TAU;

// Number of phase and relative frequency buckets
//
// This gives the spatial resolution of the final map of sawtooth wave error.
// More buckets take more time to compute, but allow studying finer details of
// the band-limited saw's error landscape.
//
const PHASE_BUCKETS: usize = 1000;
const RELATIVE_FREQ_BUCKETS: usize = 2000;

/// Minimmal number of error samples per bucket
///
/// Forcing multiple samples per bucket will produce a less noisy error map,
/// much like super-sampling antialiasing reduces aliasing artifacts in 3D
/// graphics, but the price to pay is a quadratic increase of computation time.
///
const SAMPLES_PER_BUCKET: usize = 2;

/// Size of the phase buckets
fn phase_bucket_size() -> AudioPhase {
    (PHASE_RANGE.end - PHASE_RANGE.start) / (PHASE_BUCKETS as AudioPhase)
}

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
fn log2_relative_rate_range() -> Range<AudioFrequency> {
    let start = ((SAMPLING_RATE_RANGE.start as AudioFrequency) / OSCILLATOR_FREQ_RANGE.end).log2();
    let end = ((SAMPLING_RATE_RANGE.end as AudioFrequency) / OSCILLATOR_FREQ_RANGE.start).log2();
    start..end
}

/// Size of the base-2 relative sample rate buckets
fn log2_relative_rate_bucket_size() -> AudioFrequency {
    (log2_relative_rate_range().end - log2_relative_rate_range().start)
        / (RELATIVE_FREQ_BUCKETS as AudioFrequency)
}

/// Produce a set of irregularly spaced floating-point values that span a
/// certain parameter range, such that...
///
/// - The start the parameter range if sampled.
/// - If the parameter range were divided into num_buckets regularly spaced
///   regions, at least SAMPLES_PER_BUCKET samples would fall into each bucket.
///
/// The use of irregular spacing allows us to avoid Moiré artifacts. Instead,
/// spatially unresolved high-resolution details of the sampled function will
/// manifest as random-looking noise in the sampled function's map.
///
/// The number of produced samples is random, but on average it will be twice
/// the number of buckets.
///
fn irregular_samples(range: Range<f32>, num_buckets: usize) -> impl Iterator<Item = f32> {
    genawaiter::rc::gen!({
        let mut rng = rand::thread_rng();
        let bucket_size = (range.end - range.start) / (num_buckets as f32);
        let effective_bucket_size = bucket_size / (SAMPLES_PER_BUCKET as f32);
        let mut coord = range.start;
        while coord < range.end {
            yield_!(coord);
            coord = rng.gen_range(coord..=(coord + effective_bucket_size));
        }
    })
    .into_iter()
}

/// Number of bits of error from an unlimited signal to its band-limited cousin
type ErrorBits = u8;

/// Measure how the reference saw differs from the band-unlimited saw, for a
/// certain set of parameters and at a certain phase
fn measure_error(
    sampling_rate: SamplingRateHz,
    oscillator_freq: AudioFrequency,
    phase: AudioPhase,
) -> ErrorBits {
    trace!(
        "Probing error at phase {} (= {}pi)",
        phase,
        phase / AudioPhaseMod::consts::PI
    );
    let mut oscillator = ReferenceSaw::new(sampling_rate, oscillator_freq, phase);
    let limited = oscillator.next().unwrap();
    let unlimited = unlimited_saw(phase);
    let difference = limited - unlimited;
    let correct_bits = (-(difference.abs().log2()))
        .floor()
        .max(0.0)
        .min(u8::MAX as AudioSample) as u8;
    let error_bits = (AudioSample::MANTISSA_DIGITS as u8).saturating_sub(correct_bits);
    trace!(
        "Found a difference of {} ({} bits of error)",
        difference,
        error_bits
    );
    error_bits
}

/// Map of the error of the reference saw vs the band-unlimited saw
struct ErrorMap {
    /// Map buckets in phase-major order
    error_buckets: Vec<ErrorBits>,
}
//
struct ErrorMapBucket<'error_map> {
    /// Underlying ErrorMap
    error_map: &'error_map ErrorMap,

    /// Index of the bucket within the ErrorMap
    linear_index: usize,
}
//
impl ErrorMap {
    /// Measure the error map
    fn new() -> Self {
        // Set up a bucket-filling infrastructure
        let error_buckets = std::iter::from_fn(|| Some(AtomicU8::new(0)))
            .take(RELATIVE_FREQ_BUCKETS * PHASE_BUCKETS)
            .collect::<Vec<_>>();

        // Iterate over sampling rate / oscillator frequency ratios
        let log2_relative_rates =
            irregular_samples(log2_relative_rate_range(), RELATIVE_FREQ_BUCKETS)
                .collect::<Vec<_>>();
        log2_relative_rates
            .into_par_iter()
            .for_each(|log2_relative_rate| {
                // Summon the thread-local random number generator
                let mut rng = rand::thread_rng();

                // Pick a combination of sampling rate and oscillator sampling rate
                // that matches the desired ratio.
                let relative_rate = (2.0 as AudioFrequency).powf(log2_relative_rate);
                let min_sample_rate = SAMPLING_RATE_RANGE
                    .start
                    .max((relative_rate * OSCILLATOR_FREQ_RANGE.start).ceil() as SamplingRateHz);
                let max_sample_rate = SAMPLING_RATE_RANGE
                    .end
                    .min((relative_rate * OSCILLATOR_FREQ_RANGE.end).floor() as SamplingRateHz)
                    .max(min_sample_rate);
                let sampling_rate = rng.gen_range(min_sample_rate..=max_sample_rate);
                let oscillator_freq = (sampling_rate as AudioFrequency) / relative_rate;
                debug!(
                    "Sampling a {} Hz saw at {} Hz (relative rate = {})",
                    oscillator_freq, sampling_rate, relative_rate,
                );

                // Collect data about the oscillator's error for those parameters
                for phase in irregular_samples(PHASE_RANGE, PHASE_BUCKETS) {
                    let error = measure_error(sampling_rate, oscillator_freq, phase);
                    let relative_rate_bucket = ((log2_relative_rate
                        - log2_relative_rate_range().start)
                        / log2_relative_rate_bucket_size())
                    .trunc() as usize;
                    assert!(relative_rate_bucket < RELATIVE_FREQ_BUCKETS);
                    let phase_bucket =
                        ((phase - PHASE_RANGE.start) / phase_bucket_size()).trunc() as usize;
                    assert!(phase_bucket < RELATIVE_FREQ_BUCKETS);
                    let bucket_idx = relative_rate_bucket * PHASE_BUCKETS + phase_bucket;
                    error_buckets[bucket_idx].fetch_max(error, Ordering::Relaxed);
                }
            });

        // De-atomify the buckets
        let error_buckets = error_buckets
            .into_iter()
            .map(AtomicU8::into_inner)
            .collect();
        Self { error_buckets }
    }

    /// Iterate over the buckets of the error map
    fn iter(&self) -> impl Iterator<Item = ErrorMapBucket> {
        (0..self.error_buckets.len()).map(move |linear_index| ErrorMapBucket {
            error_map: self,
            linear_index,
        })
    }
}
//
impl ErrorMapBucket<'_> {
    /// 2D bucket index
    ///
    /// Both indices are zero-based, the first index represents the sampling
    /// rate / oscillator frequency coordinate and the second index represents
    /// the phase coordinate.
    fn index(&self) -> (usize, usize) {
        let relative_rate_bucket = self.linear_index / PHASE_BUCKETS;
        let phase_bucket = self.linear_index % PHASE_BUCKETS;
        (relative_rate_bucket, phase_bucket)
    }

    /// Top-left coordinate of the bucket in the error map
    fn start(&self) -> (AudioFrequency, AudioPhase) {
        let (relative_rate_bucket, phase_bucket) = self.index();
        let log2_relative_rate = log2_relative_rate_range().start
            + (relative_rate_bucket as AudioFrequency) * log2_relative_rate_bucket_size();
        let relative_rate = (2.0 as AudioFrequency).powf(log2_relative_rate);
        let phase = PHASE_RANGE.start + (phase_bucket as AudioPhase) * phase_bucket_size();
        (relative_rate, phase)
    }

    /// Center coordinate of the bucket
    fn center(&self) -> (AudioFrequency, AudioPhase) {
        let (start_rate, start_phase) = self.start();
        let (end_rate, end_phase) = self.end();
        (
            (start_rate + end_rate) / 2.0,
            (start_phase + end_phase) / 2.0,
        )
    }

    /// Bottom-right coordinate of the bucket in the error map
    fn end(&self) -> (AudioFrequency, AudioPhase) {
        let (mut relative_rate_bucket, mut phase_bucket) = self.index();
        relative_rate_bucket += 1;
        phase_bucket += 1;
        let linear_index = relative_rate_bucket * PHASE_BUCKETS + phase_bucket;
        ErrorMapBucket {
            error_map: self.error_map,
            linear_index,
        }
        .start()
    }

    /// Measured error within that bucket
    fn error(&self) -> ErrorBits {
        self.error_map.error_buckets[self.linear_index]
    }
}

fn reference_saw_impl() -> Result<(), Box<dyn std::error::Error>> {
    use AudioPhaseMod::consts::PI;

    // Set up some infrastructure
    env_logger::init();

    // Prepare to plot the error data
    const X_MARGIN: u32 = 50;
    const Y_MARGIN: u32 = 60;
    const LABEL_SIZE: u32 = 20;
    const NUM_X_LABELS: usize = 12;
    const NUM_Y_LABELS: usize = 18;
    let relative_freq_buckets = RELATIVE_FREQ_BUCKETS as u32;
    let phase_buckets = PHASE_BUCKETS as u32;
    let phase_range_in_pi = (PHASE_RANGE.start / PI)..(PHASE_RANGE.end / PI);
    let root = BitMapBackend::new(
        "reference_saw_vs_unlimited.png",
        (relative_freq_buckets + Y_MARGIN, phase_buckets + X_MARGIN),
    )
    .into_drawing_area();
    root.fill(&WHITE)?;
    let mut chart = ChartBuilder::on(&root)
        .margin(0)
        .x_label_area_size(X_MARGIN)
        .y_label_area_size(Y_MARGIN)
        .build_cartesian_2d(log2_relative_rate_range(), phase_range_in_pi)?;
    chart
        .configure_mesh()
        .disable_mesh()
        .label_style(("sans-serif", LABEL_SIZE))
        .x_desc("Sampling rate / oscillator frequency")
        .x_label_formatter(&|&log2_x| format!("{:.0}", (2.0f32).powf(log2_x)))
        .x_labels(NUM_X_LABELS)
        .y_desc("Phase in units of pi")
        .y_labels(NUM_Y_LABELS)
        .draw()?;
    let plotting_area = chart.plotting_area();
    assert_eq!(
        plotting_area.dim_in_pixel(),
        (relative_freq_buckets, phase_buckets)
    );
    let (base_x, base_y) = plotting_area.get_base_pixel();

    // Map the error landscape
    let error_map = ErrorMap::new();

    // Draw the error map
    trace!("Error(freq, phase) table is:");
    trace!("relrate,phase,error");
    for bucket in error_map.iter() {
        // Check out the current error map bucket
        // FIXME: For this to become a real test, we should also inspect the error
        let error = bucket.error();
        let (relative_rate, phase) = bucket.center();
        trace!("{},{},{}", relative_rate, phase, error);

        // Plot the current bucket on the error map
        let scaled_error = ((error as f32 / AudioSample::MANTISSA_DIGITS as f32) * 255.0) as u8;
        let (x_idx, y_idx) = bucket.index();
        root.draw_pixel(
            (x_idx as i32 + base_x, y_idx as i32 + base_y),
            &RGBColor(scaled_error, 0, 0),
        )?;
    }
    Ok(())
}

#[test]
#[ignore]
fn reference_saw() {
    reference_saw_impl().unwrap()
}
