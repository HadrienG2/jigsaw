//! This test studies after how many iterations signals which approximate the
//! fundamental sinus become erronerous and must be reset with a fresh exact
//! sinus computation.

mod shared;

use crate::shared::parameters::{OSCILLATOR_FREQ_RANGE, PHASE_RANGE, SAMPLING_RATE_RANGE};
use jigsaw::Oscillator;
use log::trace;
use rand::Rng;

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
    let sampling_rate = 44100 /* rng.gen_range(SAMPLING_RATE_RANGE) */;
    let oscillator_freq = 2000.0 /* rng.gen_range(OSCILLATOR_FREQ_RANGE) */;
    let initial_phase = rng.gen_range(PHASE_RANGE);
    let reference = jigsaw::ReferenceSaw::new(sampling_rate, oscillator_freq, initial_phase);
    let fullit = jigsaw::FullyIterativeSaw::new(sampling_rate, oscillator_freq, initial_phase);
    println!(
        "Probing approximation breakdown point for rate={}Hz, freq={}Hz, init_phase={}pi",
        sampling_rate, oscillator_freq, initial_phase
    );
    let breakdown_iterations = reference
        .zip(fullit)
        .inspect(|(signal, reference)| trace!("Signal={}, reference={}", signal, reference))
        .position(|(signal, reference)| {
            (signal - reference).abs() > (1 << 8) as f32 * f32::EPSILON
        });
    println!(
        "Reached breakdown point after {:?} iterations",
        breakdown_iterations
    );
}
