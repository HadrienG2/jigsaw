//! Build and plot maps of some quantity across the oscillator parameter space

use crate::parameters::{
    bucket_start, irregular_samples, log2_relative_rate_range, NUM_PHASE_BUCKETS,
    NUM_RELATIVE_FREQ_BUCKETS, OSCILLATOR_FREQ_RANGE, PHASE_RANGE, SAMPLING_RATE_RANGE,
};
use core::sync::atomic::{AtomicU8, Ordering};
use jigsaw::{AudioFrequency, AudioPhase, SamplingRateHz};
use log::debug;
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
