//! Test utilities for studying errors

mod signal;

use crate::signal::Signal;
use jigsaw::{AudioFrequency, AudioPhase, AudioPhaseMod, AudioSample, SamplingRateHz};
use log::trace;

/// Number of bits of error from an unlimited signal to its band-limited cousin
pub type ErrorBits = u8;

/// Measure how the reference saw differs from the band-unlimited saw, for a
/// certain set of parameters and at a certain phase
pub fn measure_error(
    signal: &impl Signal,
    reference: &impl Signal,
    sampling_rate: SamplingRateHz,
    oscillator_freq: AudioFrequency,
    phase: AudioPhase,
) -> ErrorBits {
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
