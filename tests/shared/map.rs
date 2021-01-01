//! Build and plot maps of some quantity across the oscillator parameter space

use super::parameters::{
    bucket_start, irregular_samples, log2_relative_rate_range, NUM_PHASE_BUCKETS,
    NUM_RELATIVE_FREQ_BUCKETS, OSCILLATOR_FREQ_RANGE, PHASE_RANGE, SAMPLING_RATE_RANGE,
};
use core::sync::atomic::{AtomicU8, Ordering};
use jigsaw::{AudioFrequency, AudioPhase, AudioPhaseMod, SamplingRateHz};
use log::debug;
use plotters::prelude::*;
use rand::Rng;
use rayon::prelude::*;

/// Map of some quantity across the oscillator parameter space
pub struct OscillatorMap {
    /// Map buckets in phase-major order
    buckets: Vec<u8>,
}
//
pub struct OscillatorMapBucket<'oscillator_map> {
    /// Underlying OscillatorMap
    map: &'oscillator_map OscillatorMap,

    /// Index of the bucket within the OscillatorMap
    linear_index: usize,
}
//
impl OscillatorMap {
    /// Map a quantity across the oscillator parameter space
    pub fn measure(
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
    pub fn iter(&self) -> impl Iterator<Item = OscillatorMapBucket> {
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
    pub fn index(&self) -> (usize, usize) {
        let relative_rate_bucket = self.linear_index / NUM_PHASE_BUCKETS;
        let phase_bucket = self.linear_index % NUM_PHASE_BUCKETS;
        (relative_rate_bucket, phase_bucket)
    }

    /// Top-left coordinate of the bucket in the map
    pub fn start(&self) -> (AudioFrequency, AudioPhase) {
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
    pub fn center(&self) -> (AudioFrequency, AudioPhase) {
        let (start_rate, start_phase) = self.start();
        let (end_rate, end_phase) = self.end();
        (
            (start_rate + end_rate) / 2.0,
            (start_phase + end_phase) / 2.0,
        )
    }

    /// Bottom-right coordinate of the bucket in the map
    pub fn end(&self) -> (AudioFrequency, AudioPhase) {
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
    pub fn data(&self) -> u8 {
        self.map.buckets[self.linear_index]
    }
}

/// Given a recipe to build a map, build and plot it in such a way that failures
/// to set up the plot will be reported right away, not after the lengthy
/// process of acquiring the map's data.
///
/// Then forward the map to the caller, so that it can perform potentially
/// panicking checks after the plot has been saved to disk.
///
pub fn map_and_plot(
    map_recipe: impl FnOnce() -> OscillatorMap,
    plot_filename: &str,
    mut data_to_color: impl FnMut(u8) -> RGBColor,
) -> Result<OscillatorMap, Box<dyn std::error::Error>> {
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
    let error_map = map_recipe();

    // Draw the error map
    for bucket in error_map.iter() {
        let color = data_to_color(bucket.data());
        let (x_idx, y_idx) = bucket.index();
        root.draw_pixel((x_idx as i32 + base_x, y_idx as i32 + base_y), &color)?;
    }

    // Bubble up the map to the caller
    Ok(error_map)
}
