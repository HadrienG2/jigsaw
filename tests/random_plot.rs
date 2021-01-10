//! This integration test picks a random set of parameters and draws all
//! versions of a signal for those parameters

mod shared;

use crate::shared::{
    logger::init_logger,
    parameters::{
        irregular_samples, NUM_PHASE_BUCKETS, OSCILLATOR_FREQ_RANGE, PHASE_RANGE, SAMPLING_RATES,
    },
    signal::{BandLimitedSignal, Signal, UnlimitedSignal},
};
use rand::prelude::*;

#[test]
#[ignore]
/// Compare all saw implementations at a random point of the parameter space
fn compare_random_saws() {
    init_logger();
    let mut rng = rand::thread_rng();
    let sampling_rate = *SAMPLING_RATES.choose(&mut rng).unwrap();
    let oscillator_freq = rng.gen_range(OSCILLATOR_FREQ_RANGE);
    let unlimited = UnlimitedSignal::new(jigsaw::unlimited_saw);
    let reference = BandLimitedSignal::<jigsaw::ReferenceSaw>::new();
    let f32sin = BandLimitedSignal::<jigsaw::F32SinSaw>::new();
    let itersin = BandLimitedSignal::<jigsaw::IterativeSinSaw>::new();
    let invmul = BandLimitedSignal::<jigsaw::InvMulSaw>::new();
    let smartharms = BandLimitedSignal::<jigsaw::SmartHarmonicsSaw>::new();
    let fullit = BandLimitedSignal::<jigsaw::FullyIterativeSaw>::new();
    println!("phase,unlimited,reference,f32sin,itersin,invmul,smartharms,fullit");
    for (_phase_bucket_idx, phases) in irregular_samples(PHASE_RANGE, NUM_PHASE_BUCKETS) {
        // FIXME: Do the plot ourselves instead of printing CSV
        // FIXME: Turn this into a real test?
        for phase in phases.iter().copied() {
            println!(
                "{},{},{},{},{},{},{},{}",
                phase,
                unlimited.measure(sampling_rate, oscillator_freq, phase),
                reference.measure(sampling_rate, oscillator_freq, phase),
                f32sin.measure(sampling_rate, oscillator_freq, phase),
                itersin.measure(sampling_rate, oscillator_freq, phase),
                invmul.measure(sampling_rate, oscillator_freq, phase),
                smartharms.measure(sampling_rate, oscillator_freq, phase),
                fullit.measure(sampling_rate, oscillator_freq, phase),
            );
        }
    }
}
