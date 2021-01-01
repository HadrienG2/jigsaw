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
///
pub fn unlimited_saw(phase: AudioPhase) -> AudioSample {
    use AudioPhaseMod::consts::{PI, TAU};
    (phase + PI).rem_euclid(TAU) / PI - 1.0
}

/// Common setup for any band-limited sawtooth wave generator
///
/// Returns a suitable phase generator and the number of harmonics to be generated
///
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

/// Given an algorithm to generate the harmonics of a sinus, compute a
/// band-limited saw function with maximal precision
fn precise_saw<Harmonics, HarmonicsGenerator>(
    phase: AudioPhase,
    sin_harmonics: HarmonicsGenerator,
) -> AudioSample
where
    HarmonicsGenerator: FnOnce(f64) -> Harmonics,
    Harmonics: Iterator<Item = f64>,
{
    let phase = phase as f64 - std::f64::consts::PI;
    let mut accumulator = 0.0;
    for (idx, sin) in sin_harmonics(phase).enumerate() {
        accumulator -= sin / ((idx + 1) as f64);
    }
    accumulator /= std::f64::consts::FRAC_PI_2;
    accumulator as _
}

/// Reference sawtooth generator
///
/// This computes a band-limited sawtooth wave with maximal precision, at the
/// cost of speed. It is intended as a speed and precision reference against
/// which other speed-optimized sawtooth wave generators can be compared
///
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
        self.phase.next().map(|phase| {
            precise_saw(phase, |phase| {
                synthesis::sin_harmonics_precise::<f64>(phase, self.num_harmonics)
            })
        })
    }
}

/// Variation of ReferenceSaw that uses a single-precision sinus
///
/// This variation is about 30% faster, but it loses 10 bits of precision,
/// which means that the result would be distinguishable from that of the
/// ReferenceSaw in a 16-bit CD recording. That seems unacceptable.
///
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
        self.phase.next().map(|phase| {
            precise_saw(phase, |phase| {
                synthesis::sin_harmonics_precise::<f32>(phase, self.num_harmonics)
            })
        })
    }
}

/// Sawtooth generator that uses the sin((n+1)x) and cos((n+1)x) trigonometric
/// identities to generate all harmonics from a "seed" fundamental sine.
///
/// This algorithm produces identical results to the ReferenceSaw, but about 3x
/// faster in scenarios where its performance is bottlenecked by harmonics
/// generation (such as sampling a 20 Hz sinus at 192 kHz).
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
        self.phase.next().map(|phase| {
            precise_saw(phase, |phase| {
                synthesis::sincos_harmonics_iterative(phase.sin_cos(), self.num_harmonics)
                    .map(|(sin, _cos)| sin)
            })
        })
    }
}

/// Buffer that is suitable for storing Fourier coefficients
type FourierCoefficients = Box<[f64]>;

/// Precompute the Fourier coefficients of a band-limited saw function
fn saw_fourier_coefficients(num_harmonics: HarmonicsCounter) -> impl Iterator<Item = f64> {
    use std::f64::consts::FRAC_PI_2;
    (1..=num_harmonics).map(|harmonic| -1.0 / (harmonic as f64 * FRAC_PI_2))
}

/// Variant of IterativeSinSaw that precomputes the Fourier coefficients
///
/// For some reason, this optimization was not worthwhile on previous
/// algorithms. I suspect this to be a matter of cache locality: when the libm
/// sinus implementation is invoked, it probably flushes the cache line
/// associated with the Fourier coefficients, which could nullify the CPU
/// performance benefit of not recomputing them on every sample.
///
pub struct InvMulSaw {
    // Underlying oscillator phase iterator
    phase: OscillatorPhase,

    // Fourier coefficients of harmonics 1 and above
    fourier_coefficients: FourierCoefficients,
}
//
impl Oscillator for InvMulSaw {
    /// Set up a sawtooth oscillator.
    fn new(
        sampling_rate: SamplingRateHz,
        oscillator_freq: AudioFrequency,
        initial_phase: AudioPhase,
    ) -> Self {
        let (phase, num_harmonics) = setup_saw(sampling_rate, oscillator_freq, initial_phase);
        let fourier_coefficients = saw_fourier_coefficients(num_harmonics).collect();
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
        let mut accumulator = 0.0;
        let num_harmonics = self.fourier_coefficients.len() as HarmonicsCounter;
        for (&coeff, (sin, _cos)) in
            self.fourier_coefficients
                .iter()
                .zip(synthesis::sincos_harmonics_iterative(
                    phase.sin_cos(),
                    num_harmonics,
                ))
        {
            accumulator += coeff * sin;
        }
        Some(accumulator as _)
    }
}

/// Buffer that is suitable for storing N (sin, cos) harmonics
type SincosHarmonics = (Box<[f64]>, Box<[f64]>);

/// Build a buffer that is suitable for storing N (sin, cos) harmonics
fn sincos_harmonics_buffer(num_harmonics: HarmonicsCounter) -> SincosHarmonics {
    let sin_buf = vec![0.0; num_harmonics as usize].into_boxed_slice();
    (sin_buf.clone(), sin_buf)
}

/// Given a set of Fourier coefficients and the corresponding set of sinus
/// harmonics, synthesize an odd band-limited signal, such as a saw
fn synthesize_odd_signal(fourier_coefficients: &[f64], sinus_harmonics: &[f64]) -> AudioSample {
    // Teach the compiler that input buffers have the same size
    let num_harmonics = fourier_coefficients.len();
    let fourier_coefficients = &fourier_coefficients[..num_harmonics];
    let sinus_harmonics = &sinus_harmonics[..num_harmonics];

    // Perform signal synthesis in a manner which is amenable to automatic
    // vectorization by the compiler
    const NUM_ACCUMULATORS: usize = 1 << 6;
    let mut accumulators = [0.0; NUM_ACCUMULATORS];
    for (sins, coeffs) in sinus_harmonics
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
    accumulators[0] as _
}

/// Variant of the InvMulSaw algorithm that uses an FFT-style harmonics
/// generation algorithm for extra precision and speed.
pub struct SmartHarmonicsSaw {
    // Underlying oscillator phase iterator
    phase: OscillatorPhase,

    // Fourier coefficients of harmonics 1 and above
    fourier_coefficients: FourierCoefficients,

    // Buffers to store the sines and cosines of Fourier coefficients
    sincos_buffer: SincosHarmonics,
}
//
impl Oscillator for SmartHarmonicsSaw {
    /// Set up a sawtooth oscillator.
    fn new(
        sampling_rate: SamplingRateHz,
        oscillator_freq: AudioFrequency,
        initial_phase: AudioPhase,
    ) -> Self {
        let (phase, num_harmonics) = setup_saw(sampling_rate, oscillator_freq, initial_phase);
        let fourier_coefficients = saw_fourier_coefficients(num_harmonics).collect();
        let sincos_buffer = sincos_harmonics_buffer(num_harmonics);
        Self {
            phase,
            fourier_coefficients,
            sincos_buffer,
        }
    }
}
//
impl Iterator for SmartHarmonicsSaw {
    type Item = AudioSample;

    fn next(&mut self) -> Option<Self::Item> {
        // Teach the compiler that all our internal buffers have the same size
        let num_harmonics = self.fourier_coefficients.len();
        let fourier_coefficients = &self.fourier_coefficients[..num_harmonics];
        let sin_buf = &mut self.sincos_buffer.0[..num_harmonics];
        let cos_buf = &mut self.sincos_buffer.1[..num_harmonics];

        // Compute harmonics using an FFT-like algorithm
        let phase = self.phase.next()? as f64 - std::f64::consts::PI;
        synthesis::sincos_harmonics_smart(phase.sin_cos(), (sin_buf, cos_buf));

        // Accumulate the saw signal in a way which the compiler can vectorize
        Some(synthesize_odd_signal(fourier_coefficients, sin_buf))
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
    sincos_phase_harmonics: SincosHarmonics,

    // (sin, cos) of the underlying phase increment and its harmonics
    sincos_dphi_harmonics: SincosHarmonics,

    // Fourier coefficients of harmonics 1 and above
    fourier_coefficients: FourierCoefficients,
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
        use core::f64::consts::PI;
        let (phase, num_harmonics) = setup_saw(sampling_rate, oscillator_freq, initial_phase);

        // Compute the Fourier coefficients
        let fourier_coefficients = saw_fourier_coefficients(num_harmonics).collect();

        // Set up storage for harmonics
        let mut sincos_phase_harmonics = sincos_harmonics_buffer(num_harmonics);
        let mut sincos_dphi_harmonics = sincos_phase_harmonics.clone();
        let phase_increment = phase.phase_increment() as f64;

        // Compute the harmonics of the fundamental and the phase increment
        let fundamental_phase = initial_phase as f64 - PI - phase_increment;
        synthesis::sincos_harmonics_smart(
            fundamental_phase.sin_cos(),
            (
                &mut sincos_phase_harmonics.0[..],
                &mut sincos_phase_harmonics.1[..],
            ),
        );
        synthesis::sincos_harmonics_smart(
            phase_increment.sin_cos(),
            (
                &mut sincos_dphi_harmonics.0[..],
                &mut sincos_dphi_harmonics.1[..],
            ),
        );

        // Emit the output oscillator
        Self {
            sincos_phase_harmonics,
            sincos_dphi_harmonics,
            fourier_coefficients,
        }
    }
}
//
impl Iterator for FullyIterativeSaw {
    type Item = AudioSample;

    fn next(&mut self) -> Option<Self::Item> {
        // Teach the compiler that all our internal buffers have the same size
        let num_harmonics = self.fourier_coefficients.len();
        let fourier_coefficients = &self.fourier_coefficients[..num_harmonics];
        let sin_phase_harmonics = &mut self.sincos_phase_harmonics.0[..num_harmonics];
        let cos_phase_harmonics = &mut self.sincos_phase_harmonics.1[..num_harmonics];
        let sin_dphi_harmonics = &mut self.sincos_dphi_harmonics.0[..num_harmonics];
        let cos_dphi_harmonics = &mut self.sincos_dphi_harmonics.1[..num_harmonics];

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
        Some(synthesize_odd_signal(
            fourier_coefficients,
            sin_phase_harmonics,
        ))
    }
}

// TODO: Rework oscillators so that they accept an in-situ filter for the
//       purpose of avoiding Gibbs phenomenon when it is undesirable.
