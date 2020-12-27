mod phase;

use phase::{AudioPhase, OscillatorPhase};

/// Compute the minimal oscillator frequency that the algorithm can handle at
/// a given sampling rate.
pub fn min_oscillator_freq(sampling_rate: SamplingRateHz) -> AudioFrequency {
    assert_ne!(sampling_rate, 0, "Input sampling rate should be nonzero");
    phase::min_relative_freq() * (sampling_rate as AudioFrequency)
}

// === FOUNDATIONAL TYPES ===

/// Integer type suitable for storing an audio sampling rate in Hz
//
// In an age where high-end audio cards go to 192 000 Hz, u16 is too little, but
// u32 seems just fine for all foreseeable future.
//
pub type SamplingRateHz = u32;
#[cfg(test)]
type NonZeroSamplingRate = core::num::NonZeroU32;

/// Floating-point type suitable for storing audio frequencies
//
// Single precision enables storing frequencies with ~10^-6 precision, this is
// way better than the human ear's frequency discrimination ability.
//
pub type AudioFrequency = f32;

/// Floating-point type suitable for storing audio samples
//
// Single precision gives us a precision equivalent to 24-bit integers, which
// are currently the gold standard of integer audio sample storage, and way
// beyond human earing's amplitude discrimination abilities. This is fine.
//
pub type AudioSample = f32;

/// Floating-point type suitable for summing audio signals
//
// While AudioSample is precise enough, sums of AudioSamples can easily get
// measurably noisy. To reduce this effect, it is best to compute the sum using
// the highest available precision.
//
type AudioAccumulator = f64;

// === SAW GENERATORS ===

// TODO: Add a saw generator trait and deduplicate common functionality

// This guarantees that we can store harmonics as long as requested saw
// frequencies are above 1 Hz, which includes the whole audio range...
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
        let phase = self.phase.next()? as AudioAccumulator;
        let mut accumulator = 0.0;
        for harmonic in 1..=self.num_harmonics {
            let harmonic = harmonic as AudioAccumulator;
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
