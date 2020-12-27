//! Foundation for manipulating audio signals

use log::warn;

/// Integer type suitable for storing an audio sampling rate in Hz
//
// In an age where high-end audio cards go to 192 000 Hz, u16 is too little, but
// u32 seems just fine for all foreseeable future.
//
pub type SamplingRateHz = u32;
#[cfg(test)]
pub(crate) type NonZeroSamplingRate = core::num::NonZeroU32;

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

/// Compute the Nyquist frequency associated with an audio sampling rate
pub fn nyquist_frequency(sampling_rate: SamplingRateHz) -> AudioFrequency {
    validate_sampling_rate(sampling_rate);
    (sampling_rate as AudioFrequency) / 2.0
}

/// Validate a user-provided audio frequency
///
/// This function will panic if the input is fundamentally incorrect, and return
/// false if the input is not ideal but can be used with a reasonable chance of
/// partial success (e.g. slightly inaccurate result)..
///
pub(crate) fn validate_audio_frequency(
    (sampling_rate, freq): (SamplingRateHz, AudioFrequency),
) -> bool {
    // Check that oscillator frequency is not IEEE-754 madness.
    assert!(freq.is_finite(), "Audio frequency should be finite");

    // Warn if frequency is subnormal, as it will lead to inaccurate results
    let mut is_ideal = true;
    if !freq.is_normal() {
        is_ideal = false;
        warn!("This crate is neither designed nor tested for subnormal audio frequencies");
    }

    // Check that the requested oscillator frequency honors the
    // Shannon-Nyquist criterion.
    assert!(
        freq < (sampling_rate as AudioFrequency) / 2.0,
        "Audio frequency should honor the Shannon-Nyquist criterion"
    );
    is_ideal
}

/// Floating-point type suitable for storing audio samples
//
// Single precision gives us a precision equivalent to 24-bit integers, which
// are currently the gold standard of integer audio sample storage, and way
// beyond human earing's amplitude discrimination abilities. This is fine.
//
pub type AudioSample = f32;

#[cfg(test)]
pub(crate) mod test_tools {
    use super::*;
    use std::panic::{catch_unwind, UnwindSafe};

    /// Check that a function panics when called
    pub fn panics<R>(f: impl FnOnce() -> R + UnwindSafe) -> bool {
        catch_unwind(f).is_err()
    }

    /// Check that a certain user input passes maximally rigorous validation
    // TODO: Consider moving to a dedicated module
    pub fn is_standard<I: UnwindSafe>(
        input: I,
        validation: impl FnOnce(I) -> bool + UnwindSafe,
    ) -> bool {
        match catch_unwind(|| validation(input)) {
            Ok(true) => true,
            Ok(false) | Err(_) => false,
        }
    }

    /// We only test sampling rates above 44100 Hz because...
    /// 1. It is a minimum for perfect audio fidelity, which we should aim for
    /// 2. Few modern sound cards support less than this in hardware
    ///
    /// We do not test sampling rates which cannot be accurately stored as float
    /// because no sound card will support them for any foreseeable future.
    ///
    pub fn is_standard_rate(rate: NonZeroSamplingRate) -> bool {
        is_standard(rate.get(), validate_sampling_rate)
    }

    /// Test that a requested oscillator frequency falls into the ideal range.
    /// Other frequencies are tested via specific edge-case tests.
    pub fn is_standard_freq(rate: NonZeroSamplingRate, freq: AudioFrequency) -> bool {
        is_standard((rate.get(), freq), validate_audio_frequency)
    }
}
