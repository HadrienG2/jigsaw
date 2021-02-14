//! This test studies the error of band-limited signals at their initial phase,
//! both with respect to each other and to band-unlimited signals.

mod shared;

use crate::shared::{
    map::{map_and_plot, OscillatorMap, OscillatorMapBucket},
    signal::{BandLimitedSignal, Signal, UnlimitedSignal},
};
use core::sync::atomic::AtomicU8;
use jigsaw::{AudioSample, F64SinSaw};
use log::{info, trace};
use plotters::prelude::*;

/// Compare two implementations of a given signal and plot the difference
fn plot_initial_error(
    signal: impl Signal + Send + Sync,
    reference: impl Signal + Send + Sync,
    plot_filename: Option<&str>,
    mut error_check: impl FnMut(OscillatorMapBucket) -> bool,
) -> Result<(), Box<dyn std::error::Error>> {
    // Map the initial error, and plot it if requested to
    let make_error_map = || {
        OscillatorMap::measure(
            |sampling_rate, oscillator_freq, phase| {
                let signal = signal.measure(sampling_rate, oscillator_freq, phase);
                let reference = reference.measure(sampling_rate, oscillator_freq, phase);
                shared::error::bits_of_error(signal, reference)
            },
            AtomicU8::fetch_max,
        )
    };
    let error_map = if let Some(plot_filename) = plot_filename {
        map_and_plot(make_error_map, plot_filename, |error_bits| {
            RGBColor(
                (error_bits as u32 * 255 / AudioSample::MANTISSA_DIGITS) as u8,
                0,
                0,
            )
        })?
    } else {
        make_error_map()
    };

    // Check that the error map lives up to expectations
    trace!("Error(relfreq, phase) table is:");
    trace!("relrate,phase,error");
    let mut max_error = 0;
    for bucket in error_map.iter() {
        let error = bucket.data();
        max_error = max_error.max(error);
        let (relative_rate, phase) = bucket.center();
        trace!("{},{},{}", relative_rate, phase, error);
    }
    info!("Maximum errror is {}", max_error);
    error_map
        .iter()
        .for_each(|bucket| assert!(error_check(bucket)));

    // We're done
    Ok(())
}

#[test]
#[ignore]
/// Compare the double-precision saw to the reference saw
///
/// The two should compare equal, and since they do, we'll use the F64SinSaw as
/// our reference in subsequent tests as it is 50x faster, which means that
/// tests based on it will run a huge lot faster.
///
fn f64sin_vs_reference_saw() {
    shared::logger::init_logger();
    plot_initial_error(
        BandLimitedSignal::<F64SinSaw>::new(),
        BandLimitedSignal::<jigsaw::ReferenceSaw>::new(),
        None,
        |bucket| bucket.data() == 0,
    )
    .unwrap()
}

#[test]
#[ignore]
/// Compare the double-precision saw to a band-unlimited saw
fn f64sin_vs_unlimited_saw() {
    shared::logger::init_logger();
    plot_initial_error(
        BandLimitedSignal::<F64SinSaw>::new(),
        UnlimitedSignal::new(jigsaw::unlimited_saw),
        Some("initial_reference_vs_unlimited_saw.png"),
        |_bucket| true, // FIXME: I don't know a good error bound here
    )
    .unwrap()
}

#[test]
#[ignore]
/// Compare the saw with single-precision sinus to the double-precision saw
fn f32sin_vs_f64sin_saw() {
    shared::logger::init_logger();
    plot_initial_error(
        BandLimitedSignal::<jigsaw::F32SinSaw>::new(),
        BandLimitedSignal::<F64SinSaw>::new(),
        Some("initial_f32sin_vs_f64sin_saw.png"),
        |bucket| bucket.data() <= 10,
    )
    .unwrap()
}

#[test]
#[ignore]
/// Compare the saw with iterative sinus to the double-precision saw
fn itersin_vs_f64sin_saw() {
    shared::logger::init_logger();
    plot_initial_error(
        BandLimitedSignal::<jigsaw::IterativeSinSaw>::new(),
        BandLimitedSignal::<F64SinSaw>::new(),
        None,
        |bucket| bucket.data() == 0,
    )
    .unwrap()
}

#[test]
#[ignore]
/// Compare the saw with multiply-by-inverse to the double-precision saw
fn invmul_vs_f64sin_saw() {
    shared::logger::init_logger();
    plot_initial_error(
        BandLimitedSignal::<jigsaw::InvMulSaw>::new(),
        BandLimitedSignal::<F64SinSaw>::new(),
        None,
        |bucket| bucket.data() == 0,
    )
    .unwrap()
}

#[test]
#[ignore]
/// Compare the saw using a smart FFT-like harmonics computation method
/// to the double-precision saw
fn smartharms_vs_f64sin_saw() {
    shared::logger::init_logger();
    plot_initial_error(
        BandLimitedSignal::<jigsaw::SmartHarmonicsSaw>::new(),
        BandLimitedSignal::<F64SinSaw>::new(),
        None,
        |bucket| bucket.data() == 0,
    )
    .unwrap()
}
