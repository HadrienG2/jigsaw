pub(crate) mod audio;
mod phase;
mod synthesis;

use crate::{
    phase::{AudioPhase, OscillatorPhase},
    synthesis::HarmonicsCounter,
};

pub use crate::{
    audio::{AudioFrequency, AudioSample, SamplingRateHz, MIN_SAMPLING_RATE},
    phase::min_oscillator_freq,
    synthesis::Oscillator,
};

// === SAW GENERATORS ===

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

impl Oscillator for ReferenceSaw {
    /// Set up a sawtooth oscillator.
    fn new(
        sampling_rate: SamplingRateHz,
        oscillator_freq: AudioFrequency,
        initial_phase: AudioPhase,
    ) -> Self {
        // Validate user inputs
        audio::validate_sampling_rate(sampling_rate);
        audio::validate_audio_frequency((sampling_rate, oscillator_freq));
        // TODO: Validate initial_phase

        // Set up the phase clock
        let phase = OscillatorPhase::new(sampling_rate, oscillator_freq, initial_phase);

        // Determine how many harmonics must be generated
        let num_harmonics = synthesis::band_limited_harmonics(sampling_rate, oscillator_freq);
        synthesis::check_harmonics_precision(num_harmonics, f64::MANTISSA_DIGITS);

        // Compute the amplitude norm, which is the harmonic number
        let amplitude_norm = synthesis::harmonic_number(num_harmonics) as f64;

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
        let mut accumulator = 0.0;
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
    // TODO: Test min_x_freq functions
    #[test]
    fn it_works() {
        assert_eq!(2 + 2, 4);
    }
}
