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

/// Maximal sampling rate that this crate can accurately manipulate
//
// TODO: Constify once Rust supports it
//
pub fn max_sampling_rate() -> SamplingRateHz {
    (2 as SamplingRateHz).pow(AudioFrequency::MANTISSA_DIGITS)
}

/// Validate a user-provided sampling rate
///
/// This function will panic if the input is fundamentally incorrect, and return
/// false if the input is not ideal but can be used with a reasonable chance of
/// partial success (e.g. slightly inaccurate result)..
///
pub(crate) fn validate_sampling_rate(rate: SamplingRateHz) -> bool {
    let mut is_ideal = true;
    if rate < MIN_SAMPLING_RATE {
        is_ideal = false;
        warn!("This crate is neither designed nor tested for such low sampling rates");
    } else if rate > max_sampling_rate() {
        is_ideal = false;
        warn!("This sampling rate cannot be honored exactly and will be rounded off");
    }
    assert_ne!(rate, 0, "Input sampling rate should be nonzero");
    is_ideal
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
