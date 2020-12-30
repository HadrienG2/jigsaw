//! Mechanism for generating oscillator phase signals

use crate::{audio, AudioFrequency, SamplingRateHz};
use core::num::FpCategory;
use log::warn;
use more_asserts::*;

/// Minimum oscillator frequency that the phase generation algorithm can handle,
/// normalized by the sampling rate in use.
///
/// This limitation comes from the fact that for a precise computation, the
/// sample counter must be exactly convertible to float, but the smaller the
/// frequency is, the larget the sample counter wraparound cycle becomes...
//
// TODO: Constify this once Rust allows for it
pub fn min_relative_freq() -> AudioFrequency {
    (2.0 as AudioFrequency).powi(-(AudioPhase::MANTISSA_DIGITS as i32))
}

/// Minimum oscillator frequency that the phase generation algorithm can handle
/// at a given sampling rate.
//
// TODO: Constify this once Rust allows for it
pub fn min_oscillator_freq(sampling_rate: SamplingRateHz) -> AudioFrequency {
    audio::validate_sampling_rate(sampling_rate);
    min_relative_freq() * (sampling_rate as AudioFrequency)
}

/// Floating-point type suitable for storing an audio signal's phase
//
// Single precision enables storing audio phases with 10^-7 precision. Through
// a simple Taylor expansion, we can tell that the output harmonics
// sin(n*phi) / n are at worst precise up to dsin = n*dphi / n = dphi ~ 10^-7,
// which is fine for the same reason that storing the audio signal at 32-bit is
// fine : it gives us the best resolution that a 24-bit high end sound card is
// capable of rendrering.
//
pub type AudioPhase = f32;
pub use std::f32 as AudioPhaseMod;

/// Validate a user-provided audio phase
///
/// This function will panic if the input is fundamentally incorrect, and return
/// false if the input is not ideal but can be used with a reasonable chance of
/// partial success (e.g. slightly inaccurate result)..
///
pub(crate) fn validate_audio_phase(phase: AudioPhase) -> bool {
    match phase.classify() {
        FpCategory::Nan | FpCategory::Infinite => panic!("Oscillator phase should be finite"),
        FpCategory::Subnormal => {
            warn!("This crate is neither designed nor tested for subnormal oscillator phases");
            false
        }
        FpCategory::Normal | FpCategory::Zero => true,
    }
}

/// The sample index is represented as a floating-point number of the same type
/// as the phase to avoid conversions.
///
/// This is not a problem because for precision reason, we must anyhow operate
/// in the sample count range where sample indices are exactly representable by
/// that floating-point number type.
///
type SampleCounter = AudioPhase;

/// Simulate an oscillator's phase signal as a stream of audio-like samples
///
/// To avoid floating-point accuracy issues, the phase signal will be reset
/// after an unspecified number of 2*pi cycles and should only be relied on in
/// "modulo 2*pi" sorts of operations, like sinus computations.
///
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct OscillatorPhase {
    // Initial oscillator phase
    phase_offset: AudioPhase,

    // Phase increment on every sample
    phase_increment: AudioPhase,

    // Current phase sample
    //
    // This is an integer, but it is encoded as a floating-point to avoid
    // int -> float conversions. For this phase generator design to work out, we
    // must anyhow work in the range where the AudioPhase floating-point type
    // can exactly represent integer sample indices.
    //
    sample_idx: SampleCounter,

    // Looping point where we reset the sample index to 0 (and thus the sample
    // phase to its initial value) to avoid going out of the range where sample
    // indices are exactly representable by SampleCounter.
    //
    // This corresponds to a phase that is an exact multiple of 2*pi, so the
    // process should be largely transparent to the user.
    //
    sample_idx_cycle: SampleCounter,
}

impl OscillatorPhase {
    /// Set up a phase clock for a certain kind of oscillator
    pub fn new(
        sampling_rate: SamplingRateHz,
        oscillator_freq: AudioFrequency,
        initial_phase: AudioPhase,
    ) -> Self {
        // Check that the input parameters make sense
        audio::validate_sampling_rate(sampling_rate);
        audio::validate_audio_frequency((sampling_rate, oscillator_freq));
        validate_audio_phase(initial_phase);

        // Pick a relative operating frequency
        let relative_freq = Self::relative_freq(sampling_rate, oscillator_freq);

        // Normalize phase offset to the [0; 2*pi] range for maximal accuracy.
        let phase_offset = initial_phase % AudioPhaseMod::consts::TAU;

        // Handle the zero-frequency edge case right away
        if relative_freq == 0.0 {
            return Self {
                phase_offset,
                phase_increment: 0.0,
                sample_idx: 0.0,
                sample_idx_cycle: 1.0,
            };
        }

        // Compute the relative oscillator period, i.e. how many fractional
        // audio samples are produced during a 2*pi oscillator cycle.
        let relative_period = 1.0 / relative_freq;

        // Decode this number's IEEE-754 representation
        assert_eq!(
            AudioFrequency::RADIX,
            2,
            "This code uses base-2 IEEE-754 bit tricks"
        );
        let rel_period_ieee754 = relative_period.to_bits();
        let rel_period_fraction_bits = AudioFrequency::MANTISSA_DIGITS - 1;
        let rel_period_fraction_mask = (1 << rel_period_fraction_bits) - 1;
        let rel_period_fraction = rel_period_ieee754 & rel_period_fraction_mask;
        let rel_period_bits = std::mem::size_of::<AudioFrequency>() as u32 * 8;
        let rel_period_exponent_bits = rel_period_bits - rel_period_fraction_bits - 1;
        let rel_period_exponent_mask =
            ((1 << rel_period_exponent_bits) - 1) << rel_period_fraction_bits;
        let rel_period_biased_exponent =
            (rel_period_ieee754 & rel_period_exponent_mask) >> rel_period_fraction_bits;
        let rel_period_exponent = rel_period_biased_exponent as i32 + AudioFrequency::MIN_EXP - 2;

        // Find the number of trailing zeroes in its mantissa's fraction
        let rel_period_trailing_zeros = rel_period_fraction
            .trailing_zeros()
            .min(rel_period_fraction_bits);

        // From that and the exponent, we can deduce by which power of 2
        // relative_period must be multiplied to yield an integer, i.e. after
        // how many periods the oscillator's phase is an exact multiple of 2*pi.
        let rel_period_significant_bits = rel_period_fraction_bits - rel_period_trailing_zeros;
        let rel_periods_fractional_bits =
            (rel_period_significant_bits as i32 - rel_period_exponent).max(0);
        let rel_periods_per_cycle = (2.0 as AudioPhase).powi(rel_periods_fractional_bits as _);

        // Use that knowledge to deduce on which sample it is safe to reset the
        // oscillator's sample counter.
        let sample_idx_cycle = rel_periods_per_cycle * relative_period;
        debug_assert_eq!(
            sample_idx_cycle.fract(),
            0.0,
            "Oscillator reset period was not accurately computed, it should have been"
        );
        debug_assert_le!(
            sample_idx_cycle,
            (2.0 as AudioPhase).powi(AudioPhase::MANTISSA_DIGITS as i32 + 1),
            "Sample counter is not exactly representable by AudioPhase, it should be",
        );

        // Next, we can compute by which quantity the oscillator phase should
        // increase on every simulated sample...
        let phase_increment = AudioPhaseMod::consts::TAU * relative_freq;

        // ...and we're done !
        Self {
            phase_offset,
            phase_increment,
            sample_idx: 0.0,
            sample_idx_cycle,
        }
    }

    /// Pick a relative operating frequency
    fn relative_freq(
        sampling_rate: SamplingRateHz,
        oscillator_freq: AudioFrequency,
    ) -> AudioFrequency {
        // Check that the input parameters make sense
        audio::validate_sampling_rate(sampling_rate);
        audio::validate_audio_frequency((sampling_rate, oscillator_freq));

        // Check that the oscillator frequency can be handled by our internal
        // processing, if not adjust to the closest supported one.
        //
        // In case you're curious, the problem here is that self.sample_idx must
        // be exactly convertible to the AudioPhase floating-point type, or else
        // phases near the end of the oscillator's cycle will be inaccurate.
        // When you do the math of checking what is the breakdown point, you end
        // up on this minimal frequency criterion.
        //
        let sampling_rate = sampling_rate as AudioFrequency;
        let mut relative_freq = oscillator_freq / sampling_rate;
        let min_relative_freq = min_relative_freq();
        if relative_freq != 0.0 && relative_freq.abs() < min_relative_freq {
            warn!("Oscillator frequency cannot be honored exactly and will be rounded");
            if relative_freq.abs() >= min_relative_freq / 2.0 {
                relative_freq = relative_freq.signum() * min_relative_freq;
            } else {
                relative_freq = relative_freq.signum() * 0.0;
            }
        }
        relative_freq
    }

    /// Tell what the phase increment is
    pub fn phase_increment(
        sampling_rate: SamplingRateHz,
        oscillator_freq: AudioFrequency,
    ) -> AudioPhase {
        // Check that the input parameters make sense
        audio::validate_sampling_rate(sampling_rate);
        audio::validate_audio_frequency((sampling_rate, oscillator_freq));

        // Compute the phase increment
        AudioPhaseMod::consts::TAU * Self::relative_freq(sampling_rate, oscillator_freq)
    }
}

impl Iterator for OscillatorPhase {
    type Item = AudioPhase;

    fn next(&mut self) -> Option<Self::Item> {
        // TODO: If we manage to optimize so much that the modulo of the sample
        //       index wraparound becomes visibly expensive, provide a batched
        //       interface which generates N signal samples with a single test.
        let phase = self.phase_offset + self.sample_idx * self.phase_increment;
        self.sample_idx = (self.sample_idx + 1.0) % self.sample_idx_cycle;
        Some(phase)
    }
}

#[cfg(test)]
mod tests {
    use super::{AudioPhaseMod::consts::TAU, *};
    use crate::{
        audio::{test_tools as audio_tests, NonZeroSamplingRate},
        test_tools::{is_standard, panics},
    };
    use audio_tests::is_standard_rate;
    use quickcheck::{quickcheck, TestResult};

    /// Test that a requested oscillator frequency falls into the ideal range.
    /// Other frequencies are tested via specific edge-case tests.
    fn is_standard_freq(rate: NonZeroSamplingRate, freq: AudioFrequency) -> bool {
        audio_tests::is_standard_freq(rate, freq) && freq >= min_oscillator_freq(rate.get())
    }

    /// Test that a requested oscillator phase falls into the ideal range.
    /// Other phases are tested via specific edge-case tests.
    fn is_standard_offset(offset: AudioPhase) -> bool {
        is_standard(offset, validate_audio_phase)
    }

    /// Test that an oscillator's phase has expected initial state and behavior
    ///
    /// This test is only appropriate for configurations where the
    /// OscillatorPhase can be successfully constructed. It constructs it for
    /// you, tests all relevant properties, and returns the result.
    fn test_phase(
        rate: SamplingRateHz,
        freq: AudioFrequency,
        offset: AudioPhase,
    ) -> OscillatorPhase {
        // Make sure the right initial state was produced
        let initial_phase = OscillatorPhase::new(rate, freq, offset);
        let rate_flt = rate as AudioFrequency;
        assert_eq!(initial_phase.phase_offset, offset % TAU);
        assert_eq!(initial_phase.phase_increment, TAU * freq / rate_flt);
        assert_eq!(initial_phase.sample_idx, 0.0);
        assert_eq!(
            (initial_phase.sample_idx_cycle * initial_phase.phase_increment).fract() % TAU,
            0.0
        );

        // Make sure the output is right for this initial state
        let mut phase = initial_phase;
        let sample_idx_cycle = phase.sample_idx_cycle as u32;
        for i in 0..2 * sample_idx_cycle {
            let expected = phase.phase_offset + (i as AudioPhase) * phase.phase_increment;
            assert_eq!(phase.next(), Some(expected % TAU));
        }

        // Return the OscillatorPhase for further tests
        initial_phase
    }

    quickcheck! {
        /// Test reasonable sampling rates, frequencies, and offset values
        fn general_case(rate: NonZeroSamplingRate,
                        freq: AudioFrequency,
                        offset: AudioPhase) -> TestResult {
            // Avoid edge cases which are left to dedicated tests
            if !(is_standard_rate(rate) && is_standard_freq(rate, freq) && is_standard_offset(offset)) {
                return TestResult::discard();
            }

            // Build an oscillator and check it matches expectations
            test_phase(rate.get(), freq, offset);
            TestResult::passed()
        }

        /// Test zero sampling rate edge cases
        fn zero_sampling(freq: AudioFrequency, offset: AudioPhase) -> bool {
            panics(|| min_oscillator_freq(0))
                && panics(|| OscillatorPhase::new(0, freq, offset))
        }

        /// Test non finite oscillator frequency edge case
        fn non_finite_freq(rate: NonZeroSamplingRate, offset: AudioPhase) -> bool {
            let rate = rate.get();
            panics(|| OscillatorPhase::new(rate, AudioFrequency::NAN, offset))
                && panics(|| OscillatorPhase::new(rate, AudioFrequency::INFINITY, offset))
                && panics(|| OscillatorPhase::new(rate, AudioFrequency::NEG_INFINITY, offset))
        }

        /// Test tiny frequency edge case
        fn tiny_freqs(rate: NonZeroSamplingRate, offset: AudioPhase) -> TestResult {
            // Avoid edge cases which are left to dedicated tests
            if !(is_standard_rate(rate) && is_standard_offset(offset)) {
                return TestResult::discard();
            }

            // Build and test the oscillator
            let rate = rate.get();
            let min_oscillator_freq = min_oscillator_freq(rate);
            let small = test_phase(rate, min_oscillator_freq, offset);
            let zero = test_phase(rate, 0.0, offset);
            assert_eq!(OscillatorPhase::new(rate, 0.99 * min_oscillator_freq, offset), small);
            assert_eq!(OscillatorPhase::new(rate, 0.5 * min_oscillator_freq, offset), small);
            assert_eq!(OscillatorPhase::new(rate, 0.49 * min_oscillator_freq, offset), zero);
            assert_eq!(OscillatorPhase::new(rate, 0.01 * min_oscillator_freq, offset), zero);
            TestResult::passed()
        }

        /// Test non finite phase edge case
        fn non_finite_offset(rate: NonZeroSamplingRate, freq: AudioFrequency) -> bool {
            let rate = rate.get();
            panics(|| OscillatorPhase::new(rate, freq, AudioPhase::NAN))
                && panics(|| OscillatorPhase::new(rate, freq, AudioPhase::INFINITY))
                && panics(|| OscillatorPhase::new(rate, freq, AudioPhase::NEG_INFINITY))
        }

        /// Test subnormal offset edge case
        fn subnormal_offset(rate: NonZeroSamplingRate, freq: AudioFrequency) -> TestResult {
            // Avoid edge cases which are left to dedicated tests
            if !(is_standard_rate(rate) && is_standard_freq(rate, freq)) {
                return TestResult::discard();
            }

            // Check that subnormal offsets are handled well
            let rate = rate.get();
            test_phase(rate, freq, 0.99 * AudioPhase::MIN_POSITIVE);
            test_phase(rate, freq, 0.01 * AudioPhase::MIN_POSITIVE);
            TestResult::passed()
        }
    }
}
