//! This test studies after how many iterations signals which approximate the
//! fundamental sinus become erronerous and must be reset with a fresh exact
//! sinus computation.

mod shared;

use crate::shared::{
    error::bits_of_error,
    parameters::{OSCILLATOR_FREQ_RANGE, PHASE_RANGE, SAMPLING_RATE_RANGE},
};
use jigsaw::{AudioSample, Oscillator};
use log::{debug, info, trace};
use rand::Rng;

/// Number of iterations before a signal with cumulative error is considered
/// stable with respect to its reference
const MAX_ITERATIONS: usize = 1 << 32;

/// Number of bits of error between a signal and its reference after which we
/// don't really care how much worse it's going to get.
const MAX_ERROR_BITS: u8 = AudioSample::MANTISSA_DIGITS as u8 - 16;

// TODO: Measure and map the number of iterations that it takes before a signal
//       with cumulative error becomes different from its reference by more than
//       N bits of error.

/// Draft of a future test to study cumulative error in signals
fn study_cumulative_error(
    signal: impl Iterator<Item = AudioSample>,
    reference: impl Iterator<Item = AudioSample>,
) {
    let mut max_error_bits = 0;
    for (idx, (signal, reference)) in signal.zip(reference).take(MAX_ITERATIONS).enumerate() {
        let error_bits = bits_of_error(signal, reference);
        trace!(
            "Signal={}, reference={} ({} bits of error)",
            signal,
            reference,
            error_bits
        );
        if error_bits > max_error_bits {
            debug!(
                "First observed {} bits of error after {} iteration(s) ({} vs {})",
                error_bits, idx, signal, reference
            );
            max_error_bits = error_bits;
            if error_bits >= MAX_ERROR_BITS {
                info!("Reached maximal error after {} iteration(s)", idx);
                return;
            }
        }
    }
    info!("Reached the iteration limit, but is the signal really stable?");
}

// TODO: Generalize by...
// - Studying how this changes wrt relative frequency
// - Studying reference vs unlimited as well
// ...and please draw the fucking plot
#[test]
#[ignore]
/// Compare the iterative saw to the reference, check after how many iterations
/// it starts producing different results.
fn search_cumulative_error_breakdown() {
    shared::logger::init_logger();
    let mut rng = rand::thread_rng();
    for _ in 0..100 {
        let sampling_rate = rng.gen_range(SAMPLING_RATE_RANGE);
        let oscillator_freq = rng.gen_range(OSCILLATOR_FREQ_RANGE);
        let initial_phase = rng.gen_range(PHASE_RANGE);
        let reference = jigsaw::ReferenceSaw::new(sampling_rate, oscillator_freq, initial_phase);
        let fullit = jigsaw::FullyIterativeSaw::new(sampling_rate, oscillator_freq, initial_phase);
        info!(
            "Probing approximation breakdown point for rate={}Hz, freq={}Hz, init_phase={}pi",
            sampling_rate,
            oscillator_freq,
            initial_phase / jigsaw::AudioPhaseMod::consts::PI,
        );
        let relative_freq = oscillator_freq / (sampling_rate as f32);
        debug!(
            "rel.freq={}, harm.count={}, sincos(init.phi)={:?}, sincos(dphi)={:?}",
            relative_freq,
            ((sampling_rate as f32) / (2.0 * oscillator_freq)).trunc(),
            initial_phase.sin_cos(),
            (jigsaw::AudioPhaseMod::consts::TAU * relative_freq).sin_cos()
        );
        study_cumulative_error(fullit, reference);
    }
}
