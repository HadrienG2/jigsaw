//! This integration test picks a random set of parameters and draws all
//! versions of a signal for those parameters

mod logger;
mod parameters;
mod signal;

use crate::{
    parameters::{
        irregular_samples, NUM_PHASE_BUCKETS, OSCILLATOR_FREQ_RANGE, PHASE_RANGE,
        SAMPLING_RATE_RANGE,
    },
    signal::{BandLimitedSignal, Signal, UnlimitedSignal},
};
use rand::Rng;

#[test]
#[ignore]
/// Compare the optimized saw to the reference saw
fn compare_saws() {
    logger::init_logger();
    let mut rng = rand::thread_rng();
    let sampling_rate = rng.gen_range(SAMPLING_RATE_RANGE);
    let oscillator_freq = rng.gen_range(OSCILLATOR_FREQ_RANGE);
    let unlimited = UnlimitedSignal::new(jigsaw::unlimited_saw);
    let reference = BandLimitedSignal::<jigsaw::ReferenceSaw>::new();
    let f32sin = BandLimitedSignal::<jigsaw::F32SinSaw>::new();
    let itersin = BandLimitedSignal::<jigsaw::IterativeSinSaw>::new();
    let invmul = BandLimitedSignal::<jigsaw::InvMulSaw>::new();
    let smartharms = BandLimitedSignal::<jigsaw::SmartHarmonicsSaw>::new();
    println!("phase,unlimited,reference,f32sin,itersin,invmul,smartharms");
    for (_phase_bucket_idx, phases) in irregular_samples(PHASE_RANGE, NUM_PHASE_BUCKETS) {
        // FIXME: Do the plot ourselves instead of printing CSV
        // FIXME: Turn this into a real test
        for phase in phases.iter().copied() {
            println!(
                "{},{},{},{},{},{},{}",
                phase,
                unlimited.measure(sampling_rate, oscillator_freq, phase),
                reference.measure(sampling_rate, oscillator_freq, phase),
                f32sin.measure(sampling_rate, oscillator_freq, phase),
                itersin.measure(sampling_rate, oscillator_freq, phase),
                invmul.measure(sampling_rate, oscillator_freq, phase),
                smartharms.measure(sampling_rate, oscillator_freq, phase),
            );
        }
    }
}
