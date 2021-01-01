//! Test utilities for studying differences between signals

use crate::signal::Signal;
use jigsaw::{AudioFrequency, AudioPhase, AudioPhaseMod, AudioSample, SamplingRateHz};
use log::trace;

/// Measure how a signal differs from another, for a certain set of parameters
/// and at a certain phase
pub fn measure_error(
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
