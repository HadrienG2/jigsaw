pub(crate) mod audio;
mod phase;
mod synthesis;
#[cfg(test)]
pub(crate) mod test_tools;

use crate::synthesis::HarmonicsCounter;

pub use crate::{
    audio::{AudioFrequency, AudioSample, SamplingRateHz, MAX_SAMPLING_RATE, MIN_SAMPLING_RATE},
    phase::{min_oscillator_freq, min_relative_freq, AudioPhase, AudioPhaseMod},
    synthesis::Oscillator,
};

#[cfg(not(test))]
pub use crate::phase::OscillatorPhase;
#[cfg(test)]
pub use crate::phase::OscillatorPhase;

// === SAW GENERATORS ===

/// Sawtooth wave without band limiting
///
/// It is not correct to directly sample this ideal wave, in the sense that
/// sending a regularly sampled version of this wave through the digital-analog
/// converter of a sound card will result in spectral aliasing artifacts.
///
/// However, this waveform can be used for validation purposes, that's what the
/// band_limiting_error integration test does.
pub fn unlimited_saw(phase: AudioPhase) -> AudioSample {
    use AudioPhaseMod::consts::{PI, TAU};
    (phase + PI).rem_euclid(TAU) / PI - 1.0
}

// This computes a band-limited sawtooth wave with maximal precision, at the
// cost of speed. It is intended as a speed and precision reference against
// which other speed-optimized sawtooth wave generators can be compared
pub struct ReferenceSaw {
    // Underlying oscillator phase iterator
    phase: OscillatorPhase,

    // Number of harmonics to be generated
    num_harmonics: HarmonicsCounter,
}
//
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
//
impl Iterator for ReferenceSaw {
    type Item = AudioSample;

    fn next(&mut self) -> Option<Self::Item> {
        // This implementation is meant to be as accurate as possible, not fast.
        // So it uses double precision and avoids "performance over precision"
        // tricks such as turning division into multiplication by inverse.
        let phase = self.phase.next()? as f64 - std::f64::consts::PI;
        let mut accumulator = 0.0;
        // FIXME: This expression is wrong and must be fixed
        for harmonic in 1..=self.num_harmonics {
            let harmonic = harmonic as f64;
            accumulator -= (harmonic * phase).sin() / harmonic;
        }
        accumulator /= std::f64::consts::FRAC_PI_2;
        Some(accumulator as _)
    }
}

// This is a performance-optimized version of the ReferenceSaw
//
// TODO: Deduplicate implementation wrt ReferenceSaw
pub struct OptimizedSaw {
    // Underlying oscillator phase iterator
    phase: OscillatorPhase,

    // Amplitude of the harmonics to be generated
    harmonics_amplitude: Box<[f64]>,
}
//
impl Oscillator for OptimizedSaw {
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
        let harmonics_amplitude = (1..=num_harmonics).map(|n| 1.0 / (n as f64)).collect();

        // We're ready to generate signal
        Self {
            phase,
            harmonics_amplitude,
        }
    }
}
//
impl Iterator for OptimizedSaw {
    type Item = AudioSample;

    fn next(&mut self) -> Option<Self::Item> {
        // This implementation is meant to be as accurate as possible, not fast.
        // So it uses double precision and avoids "performance over precision"
        // tricks such as turning division into multiplication by inverse.
        let phase = self.phase.next()? as f64 - std::f64::consts::PI;
        let mut accumulator = 0.0;
        for (idx, amplitude) in self.harmonics_amplitude.iter().copied().enumerate() {
            let harmonic = (idx + 1) as f64;
            accumulator -= (harmonic * phase).sin() * amplitude;
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
