pub(crate) mod audio;
mod phase;
mod synthesis;
#[cfg(test)]
pub(crate) mod test_tools;

use crate::synthesis::HarmonicsCounter;
use log::trace;

pub use crate::{
    audio::{AudioFrequency, AudioSample, SamplingRateHz, MAX_SAMPLING_RATE, MIN_SAMPLING_RATE},
    phase::{min_oscillator_freq, min_relative_freq, AudioPhase, AudioPhaseMod},
    synthesis::Oscillator,
};

#[cfg(not(test))]
pub use crate::phase::OscillatorPhase;
#[cfg(test)]
pub use crate::phase::OscillatorPhase;

/// Display diagnostic information about wave harmonics generation
const HARMONICS_DIAGNOSTICS: bool = false;

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

/// Common setup for any band-limited sawtooth wave generator
///
/// Returns a suitable phase generator and the number of harmonics to be generated
fn setup_saw(
    sampling_rate: SamplingRateHz,
    oscillator_freq: AudioFrequency,
    initial_phase: AudioPhase,
) -> (OscillatorPhase, HarmonicsCounter) {
    // Validate user inputs
    audio::validate_sampling_rate(sampling_rate);
    audio::validate_audio_frequency((sampling_rate, oscillator_freq));
    phase::validate_audio_phase(initial_phase);

    // Set up the phase clock
    let phase = OscillatorPhase::new(sampling_rate, oscillator_freq, initial_phase);

    // Determine how many harmonics must be generated
    let num_harmonics = synthesis::band_limited_harmonics(sampling_rate, oscillator_freq);
    synthesis::check_harmonics_precision(num_harmonics, f64::MANTISSA_DIGITS);

    // Emit the basic setup
    (phase, num_harmonics)
}

/// Reference sawtooth generator
///
/// This computes a band-limited sawtooth wave with maximal precision, at the
/// cost of speed. It is intended as a speed and precision reference against
/// which other speed-optimized sawtooth wave generators can be compared
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
        let (phase, num_harmonics) = setup_saw(sampling_rate, oscillator_freq, initial_phase);
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
        let phase = self.phase.next()? as f64 - std::f64::consts::PI;
        let mut accumulator = 0.0;
        for harmonic in 1..=self.num_harmonics {
            let harmonic = harmonic as f64;
            accumulator -= (harmonic * phase).sin() / harmonic;
        }
        accumulator /= std::f64::consts::FRAC_PI_2;
        Some(accumulator as _)
    }
}

/// Variation of ReferenceSaw that uses a single-precision sinus
///
/// This variation is about 30% faster, but it loses 10 bits of precision,
/// which means that the result would be distinguishable from that of the
/// ReferenceSaw in a 16-bit CD recording. That seems unacceptable.
pub struct F32SinSaw {
    // Underlying oscillator phase iterator
    phase: OscillatorPhase,

    // Number of harmonics to be generated
    num_harmonics: HarmonicsCounter,
}
//
impl Oscillator for F32SinSaw {
    /// Set up a sawtooth oscillator.
    fn new(
        sampling_rate: SamplingRateHz,
        oscillator_freq: AudioFrequency,
        initial_phase: AudioPhase,
    ) -> Self {
        let (phase, num_harmonics) = setup_saw(sampling_rate, oscillator_freq, initial_phase);
        Self {
            phase,
            num_harmonics,
        }
    }
}
//
impl Iterator for F32SinSaw {
    type Item = AudioSample;

    fn next(&mut self) -> Option<Self::Item> {
        let phase = self.phase.next()? as f64 - std::f64::consts::PI;
        let mut accumulator = 0.0;
        for harmonic in 1..=self.num_harmonics {
            let harmonic = harmonic as f64;
            accumulator -= ((harmonic * phase) as f32).sin() as f64 / harmonic;
        }
        accumulator /= std::f64::consts::FRAC_PI_2;
        Some(accumulator as _)
    }
}

/// Sawtooth generator that uses the sin((n+1)x) and cos((n+1)x) trigonometric
/// identities to generate all harmonics from a "seed" fundamental sine.
///
/// This algorithm produces identical results to the ReferenceSaw, but about 3x
/// faster in scenarios where its performance is bottlenecked by harmonics
/// generation (such as sampling a 20 Hz sinus at 192 kHz).
///
/// Note that if need be, there is headroom for precision improvements in the
/// intermediate computations, by using a more sophisticated FFT-style harmonics
/// generation harmonics that also leverages the sin(2x) and cos(2x)
/// trigonometric identities.
///
pub struct IterativeSinSaw {
    // Underlying oscillator phase iterator
    phase: OscillatorPhase,

    // Number of harmonics to be generated
    num_harmonics: HarmonicsCounter,
}
//
impl Oscillator for IterativeSinSaw {
    /// Set up a sawtooth oscillator.
    fn new(
        sampling_rate: SamplingRateHz,
        oscillator_freq: AudioFrequency,
        initial_phase: AudioPhase,
    ) -> Self {
        let (phase, num_harmonics) = setup_saw(sampling_rate, oscillator_freq, initial_phase);
        Self {
            phase,
            num_harmonics,
        }
    }
}
//
impl Iterator for IterativeSinSaw {
    type Item = AudioSample;

    fn next(&mut self) -> Option<Self::Item> {
        let phase = self.phase.next()? as f64 - std::f64::consts::PI;
        let (sin_1, cos_1) = phase.sin_cos();
        let mut sincos_n = (sin_1, cos_1);
        let mut accumulator = -sin_1;
        for harmonic in 2..=self.num_harmonics {
            let harmonic = harmonic as f64;
            let (sin_n, cos_n) = sincos_n;
            sincos_n = (sin_n * cos_1 + cos_n * sin_1, cos_n * cos_1 - sin_n * sin_1);
            accumulator -= sincos_n.0 / harmonic;
            if HARMONICS_DIAGNOSTICS {
                let error = sin_n - (harmonic * phase).sin();
                let good_bits = (-error.abs().log2()).floor().max(0.0).min(u32::MAX as f64) as u32;
                let error_bits = f64::MANTISSA_DIGITS.saturating_sub(good_bits);
                trace!(
                    "f64 error on harmonic {} is {} ({} bits)",
                    harmonic,
                    error,
                    error_bits
                );
            }
        }
        accumulator /= std::f64::consts::FRAC_PI_2;
        Some(accumulator as _)
    }
}

/// Variant of the IterativeSinSaw algorithm that replaces the division by the
/// harmonic index with a multiplication by its inverse
pub struct InvMulSaw {
    // Underlying oscillator phase iterator
    phase: OscillatorPhase,

    // Fourier coefficients of harmonics 1 and above
    fourier_coefficients: Box<[f64]>,
}
//
impl Oscillator for InvMulSaw {
    /// Set up a sawtooth oscillator.
    fn new(
        sampling_rate: SamplingRateHz,
        oscillator_freq: AudioFrequency,
        initial_phase: AudioPhase,
    ) -> Self {
        use core::f64::consts::FRAC_PI_2;
        let (phase, num_harmonics) = setup_saw(sampling_rate, oscillator_freq, initial_phase);
        let fourier_coefficients = (1..=num_harmonics)
            .map(|harmonic| -1.0 / (harmonic as f64 * FRAC_PI_2))
            .collect();
        Self {
            phase,
            fourier_coefficients,
        }
    }
}
//
impl Iterator for InvMulSaw {
    type Item = AudioSample;

    fn next(&mut self) -> Option<Self::Item> {
        let phase = self.phase.next()? as f64 - std::f64::consts::PI;
        let (sin_1, cos_1) = phase.sin_cos();
        let mut sincos_n = (sin_1, cos_1);
        let mut accumulator = 0.0;
        for fourier_coeff in self.fourier_coefficients.iter().copied() {
            accumulator += sincos_n.0 * fourier_coeff;
            let (sin_n, cos_n) = sincos_n;
            sincos_n = (sin_n * cos_1 + cos_n * sin_1, cos_n * cos_1 - sin_n * sin_1);
        }
        Some(accumulator as _)
    }
}

/// Variant of the InvMulSaw algorithm that uses an FFT-style harmonics
/// generation algorithm on the ground that in addition to its precision
/// benefits, it might be more vectorizable
pub struct SmartHarmonicsSaw {
    // Underlying oscillator phase iterator
    phase: OscillatorPhase,

    // Fourier coefficients of harmonics 1 and above
    fourier_coefficients: Box<[f64]>,

    // Buffers to store the sines and cosines of Fourier coefficients
    sincos_buf: [Box<[f64]>; 2],
}
//
impl Oscillator for SmartHarmonicsSaw {
    /// Set up a sawtooth oscillator.
    fn new(
        sampling_rate: SamplingRateHz,
        oscillator_freq: AudioFrequency,
        initial_phase: AudioPhase,
    ) -> Self {
        use core::f64::consts::FRAC_PI_2;
        let (phase, num_harmonics) = setup_saw(sampling_rate, oscillator_freq, initial_phase);
        let fourier_coefficients = (1..=num_harmonics)
            .map(|harmonic| -1.0 / (harmonic as f64 * FRAC_PI_2))
            .collect();
        let sin_buf = vec![0.0; num_harmonics as usize].into_boxed_slice();
        let sincos_buf = [sin_buf.clone(), sin_buf];
        Self {
            phase,
            fourier_coefficients,
            sincos_buf,
        }
    }
}
//
impl Iterator for SmartHarmonicsSaw {
    type Item = AudioSample;

    fn next(&mut self) -> Option<Self::Item> {
        let phase = self.phase.next()? as f64 - std::f64::consts::PI;

        // Compute harmonics using an FFT-like algorithm
        let num_harmonics = self.fourier_coefficients.len();
        let (sin_buf, cos_buf) = self.sincos_buf.split_at_mut(1);
        let sin_buf = &mut sin_buf[0][..num_harmonics];
        let cos_buf = &mut cos_buf[0][..num_harmonics];
        let (sin_1, cos_1) = phase.sin_cos();
        sin_buf[0] = sin_1;
        cos_buf[0] = cos_1;
        let mut computed_harmonics = 1;
        while computed_harmonics < num_harmonics {
            let ref_cos = cos_buf[computed_harmonics - 1];
            let ref_sin = sin_buf[computed_harmonics - 1];
            let remaining_harmonics = num_harmonics - computed_harmonics;
            let incoming_harmonics = computed_harmonics.min(remaining_harmonics);
            for old_harmonic in 0..incoming_harmonics {
                let new_harmonic = old_harmonic + computed_harmonics;
                // Safety: This is safe because by construction, neither
                //         old_harmonic not new_harmonic can go above
                //         num_harmonics. It is necessary because right now,
                //         rustc's bound check elision is not smart enough and
                //         kills autovectorization, which is unacceptable.
                unsafe {
                    *sin_buf.get_unchecked_mut(new_harmonic) = sin_buf.get_unchecked(old_harmonic)
                        * ref_cos
                        + cos_buf.get_unchecked(old_harmonic) * ref_sin;
                    *cos_buf.get_unchecked_mut(new_harmonic) = cos_buf.get_unchecked(old_harmonic)
                        * ref_cos
                        - sin_buf.get_unchecked(old_harmonic) * ref_sin;
                }
            }
            computed_harmonics *= 2;
        }

        // Accumulate the saw signal in a way which the compiler can vectorize
        const NUM_ACCUMULATORS: usize = 1 << 6;
        let mut accumulators = [0.0; NUM_ACCUMULATORS];
        for (sins, coeffs) in sin_buf
            .chunks(NUM_ACCUMULATORS)
            .zip(self.fourier_coefficients.chunks(NUM_ACCUMULATORS))
        {
            for (accumulator, (&sin, &coeff)) in
                accumulators.iter_mut().zip(sins.iter().zip(coeffs.iter()))
            {
                *accumulator += sin * coeff;
            }
        }
        let mut num_accumulators = NUM_ACCUMULATORS;
        while num_accumulators > 1 {
            num_accumulators /= 2;
            for i in 0..num_accumulators {
                accumulators[i] += accumulators[i + num_accumulators];
            }
        }
        let result: f64 = accumulators[0];

        Some(result as _)
    }
}

/// Variant of the SmartHarmonicsSaw algorithm that avoids phase computations
/// by leveraging the fact that the phase increment is constant.
///
/// This introduces a cumulative error that will, after a number of iterations,
/// become non negligible. As part of this algorithm's validation procedure,
/// I'll measure after how many iterations this happens, and introduce an
/// automatic reset procedure that prevents it from happening.
///
pub struct FullyIterativeSaw {
    // Buffers to store the sines and cosines of harmonics
    sincos_phase_harmonics: (Box<[f64]>, Box<[f64]>),

    // (sin, cos) of the underlying phase increment and its harmonics
    sincos_dphi_harmonics: (Box<[f64]>, Box<[f64]>),

    // Fourier coefficients of harmonics 1 and above
    fourier_coefficients: Box<[f64]>,
}
//
impl Oscillator for FullyIterativeSaw {
    /// Set up a sawtooth oscillator.
    fn new(
        sampling_rate: SamplingRateHz,
        oscillator_freq: AudioFrequency,
        initial_phase: AudioPhase,
    ) -> Self {
        // Do the basic setup that is common to all saw generator algorithms
        use core::f64::consts::{FRAC_PI_2, PI};
        let (phase, num_harmonics) = setup_saw(sampling_rate, oscillator_freq, initial_phase);

        // Compute the Fourier coefficients
        let fourier_coefficients = (1..=num_harmonics)
            .map(|harmonic| -1.0 / (harmonic as f64 * FRAC_PI_2))
            .collect();

        // Set up storage for harmonics
        let sin_buf = vec![0.0; num_harmonics as usize].into_boxed_slice();
        let mut sincos_phase_harmonics = (sin_buf.clone(), sin_buf);
        let mut sincos_dphi_harmonics = sincos_phase_harmonics.clone();
        let phase_increment = phase.phase_increment() as f64;
        let fundamental_phase = initial_phase as f64 - PI - phase_increment;
        Self::compute_harmonics_approx(fundamental_phase.sin_cos(), &mut sincos_phase_harmonics);
        Self::compute_harmonics_approx(phase_increment.sin_cos(), &mut sincos_dphi_harmonics);

        // Emit the output oscillator
        Self {
            sincos_phase_harmonics,
            sincos_dphi_harmonics,
            fourier_coefficients,
        }
    }
}
//
impl FullyIterativeSaw {
    /// Given a fundamental sinus' phase, compute its harmonics with maximal precision
    fn compute_harmonics_of(phase: f64, (sin_buf, cos_buf): &mut (Box<[f64]>, Box<[f64]>)) {
        let num_harmonics = sin_buf.len();
        let sin_buf = &mut sin_buf[..num_harmonics];
        let cos_buf = &mut cos_buf[..num_harmonics];
        for (idx, (sin_phase, cos_phase)) in sin_buf.iter_mut().zip(cos_buf.iter_mut()).enumerate()
        {
            let harmonic = (idx + 1) as f64;
            let sincos_phase = (harmonic * phase).sin_cos();
            *sin_phase = sincos_phase.0;
            *cos_phase = sincos_phase.1;
        }
    }

    /// Given a fundamentel sinus, compute its harmonics in an FFT-like manner
    fn compute_harmonics_approx(
        (sin, cos): (f64, f64),
        (sin_buf, cos_buf): &mut (Box<[f64]>, Box<[f64]>),
    ) {
        let num_harmonics = sin_buf.len();
        let sin_buf = &mut sin_buf[..num_harmonics];
        let cos_buf = &mut cos_buf[..num_harmonics];
        sin_buf[0] = sin;
        cos_buf[0] = cos;
        let mut computed_harmonics = 1;
        while computed_harmonics < num_harmonics {
            let ref_cos = cos_buf[computed_harmonics - 1];
            let ref_sin = sin_buf[computed_harmonics - 1];
            let remaining_harmonics = num_harmonics - computed_harmonics;
            let incoming_harmonics = computed_harmonics.min(remaining_harmonics);
            for old_harmonic in 0..incoming_harmonics {
                let new_harmonic = old_harmonic + computed_harmonics;
                // Safety: This is safe because by construction, neither
                //         old_harmonic not new_harmonic can go above
                //         num_harmonics. It is necessary because right now,
                //         rustc's bound check elision is not smart enough and
                //         kills autovectorization, which is unacceptable.
                unsafe {
                    *sin_buf.get_unchecked_mut(new_harmonic) = sin_buf.get_unchecked(old_harmonic)
                        * ref_cos
                        + cos_buf.get_unchecked(old_harmonic) * ref_sin;
                    *cos_buf.get_unchecked_mut(new_harmonic) = cos_buf.get_unchecked(old_harmonic)
                        * ref_cos
                        - sin_buf.get_unchecked(old_harmonic) * ref_sin;
                }
            }
            computed_harmonics *= 2;
        }
    }
}
//
impl Iterator for FullyIterativeSaw {
    type Item = AudioSample;

    fn next(&mut self) -> Option<Self::Item> {
        // Move all harmonics forward by a phase increment
        //
        // TODO: Add a method for generating multiple oscillator samples that
        //       uses an FFT-like approach instead of iteratively moving the
        //       fundamental sine formward. This should be more precise, at the
        //       cost of using more memory / bandwidth.
        //
        //       This could be done by reworking the Oscillator interface so
        //       that it is parametrized by a number of samples, and emits
        //       arrays of samples instead of individual samples.
        //
        //       The number of samples could probably be restricted to be a
        //       power of 2, if need be. But I think we can do without.
        //
        //       Also, I'll end up with 2D arrays, for which ndarray might be
        //       an appropriate tool.
        //
        let num_harmonics = self.fourier_coefficients.len();
        let fourier_coefficients = &self.fourier_coefficients[..num_harmonics];
        let sin_phase_harmonics = &mut self.sincos_phase_harmonics.0[..num_harmonics];
        let cos_phase_harmonics = &mut self.sincos_phase_harmonics.1[..num_harmonics];
        let sin_dphi_harmonics = &mut self.sincos_dphi_harmonics.0[..num_harmonics];
        let cos_dphi_harmonics = &mut self.sincos_dphi_harmonics.1[..num_harmonics];
        for harmonic in 0..num_harmonics {
            let new_sincos_phase = [
                sin_phase_harmonics[harmonic] * cos_dphi_harmonics[harmonic]
                    + cos_phase_harmonics[harmonic] * sin_dphi_harmonics[harmonic],
                cos_phase_harmonics[harmonic] * cos_dphi_harmonics[harmonic]
                    - sin_phase_harmonics[harmonic] * sin_dphi_harmonics[harmonic],
            ];
            sin_phase_harmonics[harmonic] = new_sincos_phase[0];
            cos_phase_harmonics[harmonic] = new_sincos_phase[1];
        }

        // Accumulate the saw signal in a SIMD-friendly way
        const NUM_ACCUMULATORS: usize = 1 << 6;
        let mut accumulators = [0.0; NUM_ACCUMULATORS];
        for (sins, coeffs) in sin_phase_harmonics
            .chunks(NUM_ACCUMULATORS)
            .zip(fourier_coefficients.chunks(NUM_ACCUMULATORS))
        {
            for (accumulator, (&sin, &coeff)) in
                accumulators.iter_mut().zip(sins.iter().zip(coeffs.iter()))
            {
                *accumulator += sin * coeff;
            }
        }
        let mut num_accumulators = NUM_ACCUMULATORS;
        while num_accumulators > 1 {
            num_accumulators /= 2;
            for i in 0..num_accumulators {
                accumulators[i] += accumulators[i + num_accumulators];
            }
        }
        let result: f64 = accumulators[0];

        Some(result as _)
    }
}

// TODO: Rework oscillators so that they accept an in-situ filter for the
//       purpose of avoiding Gibbs phenomenon when it is undesirable.
