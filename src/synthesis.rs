//! Utilities for audio synthesis

use crate::{AudioFrequency, AudioPhase, AudioSample, SamplingRateHz};

/// Type suitable for storing oscillator harmonics
//
// This guarantees that we can store harmonics as long as requested saw
// frequencies are above 1 Hz, which includes the whole audio range...
pub(crate) type HarmonicsCounter = SamplingRateHz;

/// For oscillators whose analog counterpart have infinite frequency content,
/// such as sawtooth or square waves, tell how many harmonics must be generated
/// in order to get a maximally accurate band-limited signal.
pub(crate) fn band_limited_harmonics(
    sampling_rate: SamplingRateHz,
    oscillator_freq: AudioFrequency,
) -> HarmonicsCounter {
    // TODO: Check correctness of oscillator_freq
    let nyquist_frequency = (sampling_rate as AudioFrequency) / 2.0;
    let num_harmonics = (nyquist_frequency / oscillator_freq).trunc();
    assert!(
        num_harmonics < (2.0 as AudioFrequency).powi(32),
        "Oscillator frequency is too low, harmonics count overflows counter"
    );
    num_harmonics as _
}

/// Test that an harmonics count, such as the one generated by
/// band_limited_harmonics(), can be losslessly converted to a floating-point
/// type with a mantissa of N bits.
pub(crate) fn check_harmonics_precision(num_harmonics: HarmonicsCounter, mantissa_bits: u32) {
    let harmonic_counter_bits = std::mem::size_of::<HarmonicsCounter>() * 8;
    if (mantissa_bits as usize) < harmonic_counter_bits {
        assert!(
            num_harmonics <= (2 as HarmonicsCounter).pow(mantissa_bits),
            "Too many harmonics to accurately represent them as floating-point",
        );
    }
}

/// Compute the n-th harmonic number, that is, sum(1..=n, 1/n)
pub(crate) fn harmonic_number(num_harmonics: HarmonicsCounter) -> AudioSample {
    // TODO: Find out after how many harmonics this sum becomes wrong by more
    //       than 2^-24 and assert that this limit is not reached.
    (1..=num_harmonics)
        .map(|harmonic| 1.0 / (harmonic as f64))
        .sum::<f64>() as AudioSample
}

/// This crate is all about implementing digital oscillators for audio synthesis
pub trait Oscillator: Iterator<Item = AudioSample> {
    /// Set up an oscillator
    fn new(
        sampling_rate: SamplingRateHz,
        oscillator_freq: AudioFrequency,
        initial_phase: AudioPhase,
    ) -> Self;
}

#[cfg(test)]
mod tests {
    // TODO: Add correctness tests
}
