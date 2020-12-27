//! Mechanism for generating oscillator phase signals

use crate::{AudioFrequency, SamplingRateHz};
use core::num::FpCategory;
use log::warn;

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
use std::f32 as AudioPhaseMod;

/// The sample index is represented as a floating-point number of the same type
/// as the phase to avoid conversions.
///
/// This is not a problem because for precision reason, we must anyhow operate
/// in the sample count range where sample indices are exactly representable by
/// that floating-point number type.
///
type SampleCounter = AudioPhase;

/// Minimum relative oscillator frequency that can be handled by this algorithm
///
/// This limitation comes from the fact that for a precise computation, the
/// sample counter must be exactly convertible to float, but the smaller the
/// frequency is, the larget the sample counter wraparound cycle becomes...
///
// TODO: Constify this when possible.
pub(crate) fn min_relative_freq() -> AudioFrequency {
    (2.0 as AudioFrequency).powi(-(AudioPhase::MANTISSA_DIGITS as i32))
}

/// Simulate an oscillator's phase signal as a stream of audio-like samples
///
/// To avoid floating-point accuracy issues, the phase signal will be reset
/// after an unspecified number of 2*pi cycles and should only be relied on in
/// "modulo 2*pi" sorts of operations, like sinus computations.
///
#[derive(Debug, PartialEq)]
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
        // Check that the sampling rate is sensible
        assert_ne!(sampling_rate, 0, "Input sampling rate should be nonzero");
        if sampling_rate > (2 as SamplingRateHz).pow(AudioFrequency::MANTISSA_DIGITS) {
            warn!("Sampling rate cannot be honored exactly and will be rounded");
        }
        let sampling_rate = sampling_rate as AudioFrequency;

        // Check that oscillator frequency is not IEEE-754 madness. We tolerate
        // subnormal numbers as we have a coarser threshold on excessively low
        // frequencies that will be applied later on.
        assert!(
            oscillator_freq.is_finite(),
            "Oscillator frequency should be finite"
        );

        // Normalize the frequency by the sampling rate
        let mut relative_freq = oscillator_freq / sampling_rate;

        // Check that the requested oscillator frequency honors the
        // Shannon-Nyquist criterion. After all, the whole point of this crate
        // is to generate band-limited signals...
        assert!(
            relative_freq.abs() < 0.5,
            "Oscillator frequency should honor the Shannon-Nyquist criterion"
        );

        // Check that the oscillator frequency can be handled by our internal
        // processing, if not adjust to the closest supported one.
        //
        // In case you're curious, the problem here is that self.sample_idx must
        // be exactly convertible to the AudioPhase floating-point type, or else
        // phases near the end of the oscillator's cycle will be inaccurate.
        // When you do the math of checking what is the breakdown point, you end
        // up on this minimal frequency criterion.
        //
        let min_relative_freq = min_relative_freq();
        if relative_freq != 0.0 && relative_freq.abs() < min_relative_freq {
            warn!("Oscillator frequency cannot be honored exactly and will be rounded");
            if relative_freq.abs() >= min_relative_freq / 2.0 {
                relative_freq = relative_freq.signum() * min_relative_freq;
            } else {
                relative_freq = relative_freq.signum() * 0.0;
            }
        }

        // Check that the initial oscillator phase is not IEEE-754 madness.
        // Flush subnormal phases to zero to avoid the performance hit, and
        // normalize phase to the [0; 2*pi] range for maximal accuracy.
        let phase_offset = match initial_phase.classify() {
            FpCategory::Nan | FpCategory::Infinite => {
                panic!("Initial oscillator phase should be finite")
            }
            FpCategory::Subnormal => 0.0,
            FpCategory::Normal | FpCategory::Zero => initial_phase,
        } % AudioPhaseMod::consts::TAU;

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
            .max(rel_period_fraction_bits);

        // From that and the exponent, we can deduce by which power of 2
        // relative_period must be multiplied to yield an integer, i.e. after
        // how many periods the oscillator's phase is an exact multiple of 2*pi.
        let rel_period_significant_bits = rel_period_fraction_bits - rel_period_trailing_zeros;
        let rel_periods_fractional_bits =
            ((rel_period_significant_bits as i32) - rel_period_exponent).max(0);
        let rel_periods_per_cycle = (2.0 as AudioPhase).powi(rel_periods_fractional_bits as _);

        // Use that knowledge to deduce on which sample it is safe to reset the
        // oscillator's sample counter.
        let sample_idx_cycle = rel_periods_per_cycle * relative_period;
        debug_assert_eq!(
            sample_idx_cycle.fract(),
            0.0,
            "Oscillator reset period was not accurately computed, it should have been"
        );
        debug_assert!(
            sample_idx_cycle < (2.0 as AudioPhase).powi(AudioPhase::MANTISSA_DIGITS as _),
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
    use crate::{min_oscillator_freq, NonZeroSamplingRate};
    use quickcheck::quickcheck;
    use std::panic::{catch_unwind, UnwindSafe};

    /* sampling_rate: SamplingRateHz,
    oscillator_freq: AudioFrequency,
    initial_phase: AudioPhase, */

    /* phase_offset,
    phase_increment,
    sample_idx: 0.0,
    sample_idx_cycle,*/

    fn panics<R>(f: impl FnOnce() -> R + UnwindSafe) -> bool {
        catch_unwind(f).is_err()
    }

    fn is_standard_freq(rate: NonZeroSamplingRate, freq: AudioFrequency) -> bool {
        let rate_flt = rate.get() as AudioFrequency;
        freq.is_finite() && freq >= min_relative_freq() && freq < rate_flt / 2.0
    }

    fn is_standard_offset(offset: AudioPhase) -> bool {
        offset.is_normal()
    }

    quickcheck! {
        fn zero_sampling(freq: AudioFrequency, offset: AudioPhase) -> bool {
            panics(|| min_oscillator_freq(0))
                && panics(|| OscillatorPhase::new(0, freq, offset))
        }

        fn non_finite_freq(rate: NonZeroSamplingRate, offset: AudioPhase) -> bool {
            let rate = rate.get();
            panics(|| OscillatorPhase::new(rate, AudioFrequency::NAN, offset))
                && panics(|| OscillatorPhase::new(rate, AudioFrequency::INFINITY, offset))
                && panics(|| OscillatorPhase::new(rate, AudioFrequency::NEG_INFINITY, offset))
        }

        fn zero_freq(rate: NonZeroSamplingRate, offset: AudioPhase) -> bool {
            // Avoid edge cases which are left to dedicated tests
            if !is_standard_offset(offset) {
                return true;
            }

            // Build and test the oscillator
            let phase = OscillatorPhase::new(rate.get(), 0.0, offset);
            (phase.phase_offset == offset % TAU)
                && (phase.phase_increment == 0.0)
                && (phase.sample_idx == 0.0)
                && (phase.sample_idx_cycle == 1.0)
        }

        fn tiny_freq(rate: NonZeroSamplingRate, offset: AudioPhase) -> bool {
            // Avoid edge cases which are left to dedicated tests
            if !is_standard_offset(offset) {
                return true;
            }

            // Build and test the oscillator
            let min_oscillator_freq = min_oscillator_freq(rate.get());
            let small = OscillatorPhase::new(rate.get(), min_oscillator_freq, offset);
            let tiny = OscillatorPhase::new(rate.get(), 0.0, offset);
            (OscillatorPhase::new(rate.get(), 0.99 * min_oscillator_freq, offset) == small)
                && (OscillatorPhase::new(rate.get(), 0.5 * min_oscillator_freq, offset) == small)
                && (OscillatorPhase::new(rate.get(), 0.49 * min_oscillator_freq, offset) == tiny)
                && (OscillatorPhase::new(rate.get(), 0.01 * min_oscillator_freq, offset) == tiny)
        }

        fn non_finite_offset(rate: NonZeroSamplingRate, freq: AudioFrequency) -> bool {
            let rate = rate.get();
            panics(|| OscillatorPhase::new(rate, freq, AudioPhase::NAN))
                && panics(|| OscillatorPhase::new(rate, freq, AudioPhase::INFINITY))
                && panics(|| OscillatorPhase::new(rate, freq, AudioPhase::NEG_INFINITY))
        }

        fn subnormal_offset(rate: NonZeroSamplingRate, freq: AudioFrequency) -> bool {
            // Avoid edge cases which are left to dedicated tests
            if !is_standard_freq(rate, freq) {
                return true;
            }

            // Check that subnormal offsets are flushed to zero
            let rate = rate.get();
            let zero = OscillatorPhase::new(rate, freq, 0.0);
            OscillatorPhase::new(rate, freq, 0.99 * AudioPhase::MIN_POSITIVE) == zero
        }

        fn general_case(rate: NonZeroSamplingRate, freq: AudioFrequency, offset: AudioPhase) -> bool {
            // Avoid edge cases which are left to dedicated tests
            if !(is_standard_freq(rate, freq) && is_standard_offset(offset)) {
                return true;
            }

            // Build an oscillator and check it matches expectations
            let phase = OscillatorPhase::new(rate.into(), freq, offset);
            let rate_flt = rate.get() as AudioFrequency;
            (phase.phase_offset == offset % TAU)
                && (phase.phase_increment == TAU * freq / rate_flt)
                && (phase.sample_idx == 0.0)
                && ((phase.sample_idx_cycle * phase.phase_increment).fract() % TAU == 0.0)
        }

        // FIXME: In those cases where the constructor works, test iteration
    }

    #[test]
    fn it_works() {
        assert_eq!(2 + 2, 4);
    }
}
