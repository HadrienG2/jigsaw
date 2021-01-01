//! Test utilities for studying differences between signals

use super::{map::OscillatorMap, signal::Signal};
use jigsaw::{AudioFrequency, AudioPhase, AudioPhaseMod, AudioSample, SamplingRateHz};
use log::trace;
use std::sync::atomic::AtomicU8;

/// Bit of errors between audio samples from a signal and its reference, which
/// are for now assumed to be of magnitude 1 like all self-respecting signals
pub fn bits_of_error(signal: AudioSample, reference: AudioSample) -> u8 {
    let difference = signal - reference;
    let correct_bits = (-(difference.abs().log2()))
        .floor()
        .max(0.0)
        .min(u8::MAX as AudioSample) as u8;
    (AudioSample::MANTISSA_DIGITS as u8).saturating_sub(correct_bits)
}
