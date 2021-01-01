//! This test studies the error of band-limited signals at their initial phase,
//! both with respect to each other and to band-unlimited signals.

mod error;
mod logger;
mod parameters;
pub mod signal;

use crate::{
    error::measure_error,
    logger::init_logger,
    parameters::{
        bucket_start, irregular_samples, log2_relative_rate_range, NUM_PHASE_BUCKETS,
        NUM_RELATIVE_FREQ_BUCKETS, OSCILLATOR_FREQ_RANGE, PHASE_RANGE, SAMPLING_RATE_RANGE,
    },
    signal::{BandLimitedSignal, Signal, UnlimitedSignal},
};
use core::sync::atomic::{AtomicU8, Ordering};
use jigsaw::{
    AudioFrequency, AudioPhase, AudioPhaseMod, AudioSample, ReferenceSaw, SamplingRateHz,
};
use log::{debug, info, trace};
use plotters::prelude::*;
use rand::Rng;
use rayon::prelude::*;

/// Map of some quantity across the oscillator parameter space
struct OscillatorMap {
    /// Map buckets in phase-major order
    buckets: Vec<u8>,
}
//
struct OscillatorMapBucket<'oscillator_map> {
    /// Underlying OscillatorMap
    map: &'oscillator_map OscillatorMap,

    /// Index of the bucket within the OscillatorMap
    linear_index: usize,
}
//
impl OscillatorMap {
    /// Map a quantity across the oscillator parameter space
    fn measure(
        map: impl Fn(SamplingRateHz, AudioFrequency, AudioPhase) -> u8 + Send + Sync,
        reduce: impl Fn(&AtomicU8, u8, Ordering) -> u8 + Send + Sync,
    ) -> Self {
        // Set up a bucket-filling infrastructure
        let buckets = std::iter::from_fn(|| Some(AtomicU8::new(0)))
            .take(NUM_RELATIVE_FREQ_BUCKETS * NUM_PHASE_BUCKETS)
            .collect::<Vec<_>>();

        // Iterate over sampling rate / oscillator frequency ratios
        let log2_relative_rates =
            irregular_samples(log2_relative_rate_range(), NUM_RELATIVE_FREQ_BUCKETS)
                .collect::<Vec<_>>();
        log2_relative_rates.into_par_iter().for_each(
            |(relative_rate_bucket, log2_relative_rates)| {
                for log2_relative_rate in log2_relative_rates.iter().copied() {
                    // Summon the thread-local random number generator
                    let mut rng = rand::thread_rng();

                    // Pick a combination of sampling rate and oscillator sampling rate
                    // that matches the desired ratio.
                    let relative_rate = (2.0 as AudioFrequency).powf(log2_relative_rate);
                    let min_sample_rate =
                        SAMPLING_RATE_RANGE
                            .start
                            .max((relative_rate * OSCILLATOR_FREQ_RANGE.start).ceil()
                                as SamplingRateHz);
                    let max_sample_rate = SAMPLING_RATE_RANGE
                        .end
                        .min((relative_rate * OSCILLATOR_FREQ_RANGE.end).floor() as SamplingRateHz)
                        .max(min_sample_rate);
                    let sampling_rate = rng.gen_range(min_sample_rate..=max_sample_rate);
                    let oscillator_freq = (sampling_rate as AudioFrequency) / relative_rate;
                    debug!(
                        "Sampling a {} Hz signal at {} Hz (relative rate = {})",
                        oscillator_freq, sampling_rate, relative_rate,
                    );

                    // Collect data about the oscillator for those parameters
                    for (phase_bucket, phases) in irregular_samples(PHASE_RANGE, NUM_PHASE_BUCKETS)
                    {
                        for phase in phases.iter().copied() {
                            let data = map(sampling_rate, oscillator_freq, phase);
                            let bucket_idx =
                                relative_rate_bucket * NUM_PHASE_BUCKETS + phase_bucket;
                            reduce(&buckets[bucket_idx], data, Ordering::Relaxed);
                        }
                    }
                }
            },
        );

        // De-atomify the buckets
        let buckets = buckets.into_iter().map(AtomicU8::into_inner).collect();
        Self { buckets }
    }

    /// Iterate over the buckets of the map
    fn iter(&self) -> impl Iterator<Item = OscillatorMapBucket> {
        (0..self.buckets.len()).map(move |linear_index| OscillatorMapBucket {
            map: self,
            linear_index,
        })
    }
}
//
impl OscillatorMapBucket<'_> {
    /// 2D bucket index
    ///
    /// Both indices are zero-based, the first index represents the sampling
    /// rate / oscillator frequency coordinate and the second index represents
    /// the phase coordinate.
    fn index(&self) -> (usize, usize) {
        let relative_rate_bucket = self.linear_index / NUM_PHASE_BUCKETS;
        let phase_bucket = self.linear_index % NUM_PHASE_BUCKETS;
        (relative_rate_bucket, phase_bucket)
    }

    /// Top-left coordinate of the bucket in the map
    fn start(&self) -> (AudioFrequency, AudioPhase) {
        let (relative_rate_bucket, phase_bucket) = self.index();
        let log2_relative_rate = bucket_start(
            log2_relative_rate_range(),
            NUM_PHASE_BUCKETS,
            relative_rate_bucket,
        );
        let relative_rate = (2.0 as AudioFrequency).powf(log2_relative_rate);
        let phase = bucket_start(PHASE_RANGE, NUM_PHASE_BUCKETS, phase_bucket);
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

    /// Bottom-right coordinate of the bucket in the map
    fn end(&self) -> (AudioFrequency, AudioPhase) {
        let (mut relative_rate_bucket, mut phase_bucket) = self.index();
        relative_rate_bucket += 1;
        phase_bucket += 1;
        let linear_index = relative_rate_bucket * NUM_PHASE_BUCKETS + phase_bucket;
        OscillatorMapBucket {
            map: self.map,
            linear_index,
        }
        .start()
    }

    /// Measured data within that bucket
    fn data(&self) -> u8 {
        self.map.buckets[self.linear_index]
    }
}

/// Compare two implementations of a given signal
fn plot_initial_error(
    signal: impl Signal + Send + Sync,
    reference: impl Signal + Send + Sync,
    plot_filename: &str,
) -> Result<(), Box<dyn std::error::Error>> {
    // Prepare to plot the error data
    use AudioPhaseMod::consts::PI;
    const X_MARGIN: u32 = 50;
    const Y_MARGIN: u32 = 60;
    const LABEL_SIZE: u32 = 20;
    const NUM_X_LABELS: usize = 12;
    const NUM_Y_LABELS: usize = 18;
    let relative_freq_buckets = NUM_RELATIVE_FREQ_BUCKETS as u32;
    let phase_buckets = NUM_PHASE_BUCKETS as u32;
    let phase_range_in_pi = (PHASE_RANGE.start / PI)..(PHASE_RANGE.end / PI);
    let root = BitMapBackend::new(
        plot_filename,
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
    let error_map = OscillatorMap::measure(
        |sampling_rate, oscillator_freq, initial_phase| {
            measure_error(
                &signal,
                &reference,
                sampling_rate,
                oscillator_freq,
                initial_phase,
            )
        },
        AtomicU8::fetch_max,
    );

    // Draw the error map
    trace!("Error(relfreq, phase) table is:");
    trace!("relrate,phase,error");
    let mut max_error = 0;
    for bucket in error_map.iter() {
        // Check out the current error map bucket
        // FIXME: For this to become a real test, we should also inspect the error
        let error = bucket.data();
        max_error = max_error.max(error);
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
    info!("Maximum error is {} bits", max_error);
    Ok(())
}

#[test]
#[ignore]
/// Compare the reference saw to a band-unlimited saw
fn reference_vs_unlimited_saw() {
    init_logger();
    plot_initial_error(
        BandLimitedSignal::<ReferenceSaw>::new(),
        UnlimitedSignal::new(jigsaw::unlimited_saw),
        "initial_reference_vs_unlimited_saw.png",
    )
    .unwrap()
}

#[test]
#[ignore]
/// Compare the saw with single-precision sinus to the reference saw
fn f32sin_vs_reference_saw() {
    init_logger();
    plot_initial_error(
        BandLimitedSignal::<jigsaw::F32SinSaw>::new(),
        BandLimitedSignal::<ReferenceSaw>::new(),
        "initial_f32sin_vs_reference_saw.png",
    )
    .unwrap()
}

#[test]
#[ignore]
/// Compare the saw with iterative sinus to the reference saw
fn itersin_vs_reference_saw() {
    init_logger();
    plot_initial_error(
        BandLimitedSignal::<jigsaw::IterativeSinSaw>::new(),
        BandLimitedSignal::<ReferenceSaw>::new(),
        "initial_itersin_vs_reference_saw.png",
    )
    .unwrap()
}

#[test]
#[ignore]
/// Compare the saw with multiply-by-inverse to the reference saw
fn invmul_vs_reference_saw() {
    init_logger();
    plot_initial_error(
        BandLimitedSignal::<jigsaw::InvMulSaw>::new(),
        BandLimitedSignal::<ReferenceSaw>::new(),
        "initial_invmul_vs_reference_saw.png",
    )
    .unwrap()
}

#[test]
#[ignore]
/// Compare the saw using a smart FFT-like harmonics computation method
/// to the reference saw
fn smartharms_vs_reference_saw() {
    init_logger();
    plot_initial_error(
        BandLimitedSignal::<jigsaw::SmartHarmonicsSaw>::new(),
        BandLimitedSignal::<ReferenceSaw>::new(),
        "initial_smartharms_vs_reference_saw.png",
    )
    .unwrap()
}

#[test]
#[ignore]
/// Compare the saw using a fully iterative harmonics computation method to the
/// reference saw
fn fullit_vs_reference_saw() {
    init_logger();
    plot_initial_error(
        BandLimitedSignal::<jigsaw::FullyIterativeSaw>::new(),
        BandLimitedSignal::<ReferenceSaw>::new(),
        "initial_fullit_vs_reference_saw.png",
    )
    .unwrap()
}
