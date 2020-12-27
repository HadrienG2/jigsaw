//! Foundation for manipulating audio signals

use log::warn;

/// Integer type suitable for storing an audio sampling rate in Hz
//
// In an age where high-end audio cards go to 192 000 Hz, u16 is too little, but
// u32 seems just fine for all foreseeable future.
//
pub type SamplingRateHz = u32;
pub type NonZeroSamplingRate = core::num::NonZeroU32;

/// Minimal sampling rate that this crate is tested against
//
// This limit is introduced because writing code that remains numerically stable
// at very low sampling rates is harder, and the usefulness of doing this is
// debatable when hardware support for sampling rates lower than 44.1 kHz is
// very sparse in modern audio chips.
//
pub const MIN_SAMPLING_RATE: SamplingRateHz = 44100;

/// Validate a user-provided sampling rate
pub(crate) fn validate_sampling_rate(rate: SamplingRateHz) -> NonZeroSamplingRate {
    if rate < MIN_SAMPLING_RATE {
        warn!("This crate is not tested for such low sampling rates and may break");
    }
    if rate > (2 as SamplingRateHz).pow(AudioFrequency::MANTISSA_DIGITS) {
        warn!("Sampling rate cannot be honored exactly and will be rounded");
    }
    NonZeroSamplingRate::new(rate).expect("Input sampling rate should be nonzero")
}

/// Floating-point type suitable for storing audio frequencies
//
// Single precision enables storing frequencies with ~10^-7 precision, this is
// better than the human ear's frequency discrimination ability.
//
pub type AudioFrequency = f32;

/// Floating-point type suitable for storing audio samples
//
// Single precision gives us a precision equivalent to 24-bit integers, which
// are currently the gold standard of integer audio sample storage, and way
// beyond human earing's amplitude discrimination abilities. This is fine.
//
pub type AudioSample = f32;
