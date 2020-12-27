//! Utilities for audio synthesis

use crate::{audio, AudioFrequency, AudioPhase, AudioSample, SamplingRateHz};

/// Type suitable for storing oscillator harmonics
//
// This guarantees that we can store harmonics as long as requested saw
// frequencies are above 1 Hz, which includes the whole audio range...
pub(crate) type HarmonicsCounter = SamplingRateHz;

/// Number of bits of the HarmonicsCounter integer type
const HARMONICS_COUNTER_BITS: u32 = (std::mem::size_of::<HarmonicsCounter>() as u32) * 8;

/// For oscillators whose analog counterpart have infinite frequency content,
/// such as sawtooth or square waves, tell how many harmonics must be generated
/// in order to get a maximally accurate band-limited signal.
pub(crate) fn band_limited_harmonics(
    sampling_rate: SamplingRateHz,
    oscillator_freq: AudioFrequency,
) -> HarmonicsCounter {
    // Check that the inputs make sense
    audio::validate_sampling_rate(sampling_rate);
    audio::validate_audio_frequency((sampling_rate, oscillator_freq));

    // Compute the number of harmonics of a band-limited signal
    let nyquist_frequency = audio::nyquist_frequency(sampling_rate);
    let num_harmonics = (nyquist_frequency / oscillator_freq).trunc();
    assert!(
        num_harmonics <= HarmonicsCounter::MAX as AudioFrequency,
        "Oscillator frequency is too low, harmonics count overflows counter"
    );
    num_harmonics as _
}

/// Test that an harmonics count, such as the one generated by
/// band_limited_harmonics(), can be losslessly converted to a floating-point
/// type with a mantissa of N bits.
pub(crate) fn check_harmonics_precision(num_harmonics: HarmonicsCounter, mantissa_bits: u32) {
    if mantissa_bits < HARMONICS_COUNTER_BITS {
        assert!(
            num_harmonics <= (2 as HarmonicsCounter).pow(mantissa_bits),
            "Too many harmonics to accurately represent them as floating-point",
        );
    }
}

/// Compute the n-th harmonic number, that is, sum(1..=n, 1/n)
pub(crate) fn harmonic_number(num_harmonics: HarmonicsCounter) -> AudioSample {
    let prev_number = (1..num_harmonics)
        .map(|harmonic| 1.0 / (harmonic as f64))
        .sum::<f64>();
    let curr_contrib = 1.0 / (num_harmonics as f64);
    let result = prev_number + curr_contrib;
    let error = (result - prev_number) - curr_contrib;
    assert!(
        error.abs() < AudioSample::EPSILON.into(),
        "Too many harmonics to accurately compute the harmonic number"
    );
    result as AudioSample
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
    use super::*;
    use crate::audio::{test_tools as audio_tests, NonZeroSamplingRate};
    use audio_tests::{is_standard_rate, panics};

    mod band_limited_harmonics {
        use super::*;
        use quickcheck::{quickcheck, TestResult};

        /// Minimal input frequency for this function
        fn min_freq(rate: SamplingRateHz) -> AudioFrequency {
            let max_harmonics = HarmonicsCounter::MAX as AudioFrequency;
            audio::nyquist_frequency(rate) / max_harmonics
        }

        /// Test that a requested oscillator frequency falls into the ideal range.
        /// Other frequencies are tested via specific edge-case tests.
        fn is_standard_freq(rate: NonZeroSamplingRate, freq: AudioFrequency) -> bool {
            audio_tests::is_standard_freq(rate, freq) && freq <= min_freq(rate.get())
        }

        quickcheck! {
            /// Test reasonable sampling rates and frequencies
            fn general_case(rate: NonZeroSamplingRate, freq: AudioFrequency) -> TestResult {
                // Avoid edge cases which are left to dedicated tests
                if !(is_standard_rate(rate) && is_standard_freq(rate, freq)) {
                    return TestResult::discard();
                }

                // Check the amount of band-limited harmonics
                // TODO: Handle "too many harmonics scenario as a separate test
                let rate = rate.get();
                let num_harmonics = band_limited_harmonics(rate, freq);
                let num_harmonics = num_harmonics as AudioFrequency;
                let nyquist_frequency = audio::nyquist_frequency(rate);
                assert!(num_harmonics * freq <= nyquist_frequency);
                assert!((num_harmonics + 1.0) * freq > nyquist_frequency);
                TestResult::passed()
            }

            /// Test excessively low frequencies
            fn too_many_harmonics(rate: NonZeroSamplingRate) -> bool {
                let rate = rate.get();
                panics(|| band_limited_harmonics(rate, 0.99 * min_freq(rate)))
                    && panics(|| band_limited_harmonics(rate, 0.01 * min_freq(rate)))
            }
        }
    }

    // TODO: Test other functionality

    // TODO: Add correctness tests, including a generic oscillator test that
    //       takes the band-unlimited version as input and compares against it
    //       for a very high sampling rate / frequency ratio (say,
    //       2^f32::MANTISSA_BITS), where the two versions should be ~identical.
}