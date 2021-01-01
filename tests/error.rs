//! Test utilities for studying differences between signals

use crate::{map::OscillatorMap, signal::Signal};
use jigsaw::{AudioFrequency, AudioPhase, AudioPhaseMod, AudioSample, SamplingRateHz};
use log::trace;
use std::sync::atomic::AtomicU8;

/// Measure how a signal differs from another, for a certain set of parameters
/// and at a certain phase
pub fn measure_initial_error(
    signal: &impl Signal,
    reference: &impl Signal,
    sampling_rate: SamplingRateHz,
    oscillator_freq: AudioFrequency,
    phase: AudioPhase,
) -> u8 {
    trace!(
        "Probing error at phase {} (= {}pi)",
        phase,
        phase / AudioPhaseMod::consts::PI
    );
    let signal = signal.measure(sampling_rate, oscillator_freq, phase);
    let reference = reference.measure(sampling_rate, oscillator_freq, phase);
    let difference = signal - reference;
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

/// Map the result of measure_initial_error across the parameter space
pub fn map_initial_error(
    signal: impl Signal + Send + Sync,
    reference: impl Signal + Send + Sync,
) -> OscillatorMap {
    OscillatorMap::measure(
        |sampling_rate, oscillator_freq, initial_phase| {
            measure_initial_error(
                &signal,
                &reference,
                sampling_rate,
                oscillator_freq,
                initial_phase,
            )
        },
        AtomicU8::fetch_max,
    )
}
