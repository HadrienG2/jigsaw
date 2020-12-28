pub(crate) mod audio;
mod phase;
mod synthesis;
#[cfg(test)]
pub(crate) mod test_tools;

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
        phase::validate_audio_phase(initial_phase);

        // Set up the phase clock
        let phase = OscillatorPhase::new(sampling_rate, oscillator_freq, initial_phase);

        // Determine how many harmonics must be generated
        let num_harmonics = synthesis::band_limited_harmonics(sampling_rate, oscillator_freq);
        synthesis::check_harmonics_precision(num_harmonics, f64::MANTISSA_DIGITS);

        // We're ready to generate signal
        Self {
            phase,
            num_harmonics,
        }
    }
}

impl Iterator for ReferenceSaw {
    type Item = AudioSample;

    fn next(&mut self) -> Option<Self::Item> {
        // This implementation is meant to be as accurate as possible, not fast.
        // So it uses double precision and avoids "performance over precision"
        // tricks such as turning division into multiplication by inverse.
        let phase = self.phase.next()? as f64 - std::f64::consts::PI;
        let mut accumulator = 0.0;
        for harmonic in -(self.num_harmonics as i32)..0 {
            let harmonic = harmonic as f64;
            accumulator -= (harmonic * phase).sin() / harmonic;
        }
        accumulator /= std::f64::consts::FRAC_PI_2;
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

// TODO: Rework oscillators so that they accept an in-situ filter for the
//       purpose of avoiding Gibbs phenomenon when it is undesirable.

// TODO: Add correctness tests of the reference saw
// TODO: Test correctness of other saws by comparing them to the reference
#[cfg(test)]
mod tests {
    use super::*;
    use phase::AudioPhaseMod;
    use synthesis::test_tools::test_oscillator;

    /// Continuous saw signal without band limiting
    fn unlimited_saw(phase: AudioPhase) -> AudioSample {
        use AudioPhaseMod::consts::{PI, TAU};
        (phase + PI).rem_euclid(TAU) / PI - 1.0
    }

    /// Test that our reference band-limited saw matches its continuous cousin
    /// that did not receive any band limiting.
    #[test]
    fn reference_saw() {
        test_oscillator::<ReferenceSaw>(unlimited_saw);
    }
}
