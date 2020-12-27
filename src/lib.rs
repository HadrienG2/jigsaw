use core::num::FpCategory;
use log::warn;

// === FOUNDATIONAL TYPES ===

/// Floating-point type suitable for storing audio frequencies
//
// Single precision enables storing frequencies with ~10^-6 precision, this is
// way better than the human ear's frequency discrimination ability.
//
pub type AudioFrequency = f32;

/// Floating-point type suitable for storing audio phases
//
// Assuming that audio phase is wrapped around every 2*pi, this enables storing
// phase with 10^-6 precision. In the Real World, we'll need to wrap phases
// around a bit less often, but nonetheless this seems fine.
//
pub type AudioPhase = f32;
use std::f32 as AudioPhaseMod;

/// Floating-point type suitable for storing audio samples
//
// Single precision gives us a precision equivalent to 24-bit integers, which
// are currently the gold standard of integer audio sample storage, and way
// beyond human earing's amplitude discrimination abilities. This is fine.
//
pub type AudioSample = f32;

// === OSCILLATOR PHASE CLOCK ===

// TODO: This is now mature enough, move it to a module

/// Integer type suitable for storing an audio sampling rate in Hz
//
// In an age where high-end audio cards go to 192 000 Hz, u16 is too little, but
// u32 seems just fine for all foreseeable future.
//
pub type SamplingRateHz = u32;

/// The sample index is represented as a floating-point number of the same type
/// as the phase to avoid conversions.
///
/// This is not a problem because for precision reason, we must anyhow operate
/// in the sample count range where sample indices are exactly representable by
/// that floating-point number type.
///
pub type SampleCounter = AudioPhase;

/// Simulate an oscillator's phase signal as a stream of audio-like samples
///
/// To avoid floating-point accuracy issues, the phase signal will be reset
/// after an unspecified number of 2*pi cycles and should only be relied on in
/// "modulo 2*pi" sorts of operations, like sinus computations.
///
#[derive(Debug)]
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
        let min_relative_freq = (2.0 as AudioFrequency).powi(-(AudioPhase::MANTISSA_DIGITS as i32));
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
        let sample_idx_cycle = sample_idx_cycle as _;

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

// === SAW GENERATORS ===

// TODO: Add a saw generator trait and deduplicate common functionality

// This guarantees that we can store harmonics as long as requested saw
// frequencies are above 1 Hz, which includs the whole audio range...
type HarmonicsCounter = SamplingRateHz;

// --- REFERENCE SAW ---

// This computes a band-limited sawtooth wave with maximal precision, at the
// cost of speed. It is intended as a speed and precision reference against
// which other speed-optimized sawtooth wave generators can be compared
pub struct ReferenceSaw {
    // Underlying oscillator phase iterator
    phase: OscillatorPhase,

    // Number of harmonics to be generated
    num_harmonics: HarmonicsCounter,

    // Amplitude norm, chosen so that the output is between -1.0 and +1.0
    amplitude_norm: f64,
}

impl ReferenceSaw {
    /// Set up a sawtooth oscillator.
    pub fn new(
        sampling_rate: SamplingRateHz,
        oscillator_freq: AudioFrequency,
        initial_phase: AudioPhase,
    ) -> Self {
        // Set up the oscillator's phase clock
        let phase = OscillatorPhase::new(sampling_rate, oscillator_freq, initial_phase);

        // We will generate band-limited saws via harmonic decomposition of the
        // saw function. For this process to be accurate, harmonic numbers must
        // be exactly convertible to our choice of floating-point type.
        //
        // TODO: Deduplicate this check so that all saw generators can use it
        //
        let nyquist_frequency = (sampling_rate as AudioFrequency) / 2.0;
        let num_harmonics = (nyquist_frequency / oscillator_freq).trunc();
        assert!(
            num_harmonics < (2.0 as AudioPhase).powi(AudioPhase::MANTISSA_DIGITS as i32),
            "Harmonic indices must be exactly convertible to float"
        );
        let num_harmonics = num_harmonics as HarmonicsCounter;

        // Compute the amplitude norm, which is the harmonic number
        //
        // TODO: Also make this available to other saw generators
        //
        let amplitude_norm = (1..=num_harmonics)
            .map(|harmonic| 1.0 / (harmonic as f64))
            .sum();

        // We're ready to generate signal
        Self {
            phase,
            num_harmonics,
            amplitude_norm,
        }
    }
}

impl Iterator for ReferenceSaw {
    type Item = AudioSample;

    fn next(&mut self) -> Option<Self::Item> {
        // This implementation is meant to be as accurate as possible, not fast.
        // So it uses double precision and avoids "performance over precision"
        // tricks such as turning division into multiplication by inverse.
        let phase = self.phase.next()? as f64;
        let mut accumulator = 0.0f64;
        for harmonic in 1..=self.num_harmonics {
            let harmonic = harmonic as f64;
            accumulator += (harmonic * phase).sin() / harmonic;
        }
        accumulator /= self.amplitude_norm;
        Some(accumulator as _)
    }
}

// TODO: Do a version that still uses libm sinus, but uses f32 precision,
//       multiplication by inverse, and caches all the 1/harmonic inverses in a
//       lookup table. Compare precision and performance.
//
//       This should be done after extracting as many common blocks as possible
//       from the reference implementation. We're almost there!

// TODO: Implement the sinus avoidance optimizations that this project has
//       always been meant to test, compare precision and performance.

// TODO: Add correctness tests
#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
        assert_eq!(2 + 2, 4);
    }
}
