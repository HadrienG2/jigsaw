//! This test studies the error of band-limited signals at their initial phase,
//! both with respect to each other and to band-unlimited signals.

mod shared;

use crate::shared::{
    map::{map_and_plot, OscillatorMapBucket},
    signal::{BandLimitedSignal, Signal, UnlimitedSignal},
};
use jigsaw::{AudioSample, ReferenceSaw};
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
    let make_error_map = || shared::error::map_initial_error(signal, reference);
    let error_map = if let Some(plot_filename) = plot_filename {
        map_and_plot(make_error_map, plot_filename, |error| {
            let scaled_error = ((error as f32 / AudioSample::MANTISSA_DIGITS as f32) * 255.0) as u8;
            RGBColor(scaled_error, 0, 0)
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
/// Compare the reference saw to a band-unlimited saw
fn reference_vs_unlimited_saw() {
    shared::logger::init_logger();
    plot_initial_error(
        BandLimitedSignal::<ReferenceSaw>::new(),
        UnlimitedSignal::new(jigsaw::unlimited_saw),
        Some("initial_reference_vs_unlimited_saw.png"),
        |_bucket| true, // FIXME: I don't know a good error bound here
    )
    .unwrap()
}

#[test]
#[ignore]
/// Compare the saw with single-precision sinus to the reference saw
fn f32sin_vs_reference_saw() {
    shared::logger::init_logger();
    plot_initial_error(
        BandLimitedSignal::<jigsaw::F32SinSaw>::new(),
        BandLimitedSignal::<ReferenceSaw>::new(),
        Some("initial_f32sin_vs_reference_saw.png"),
        |bucket| bucket.data() <= 10,
    )
    .unwrap()
}

#[test]
#[ignore]
/// Compare the saw with iterative sinus to the reference saw
fn itersin_vs_reference_saw() {
    shared::logger::init_logger();
    plot_initial_error(
        BandLimitedSignal::<jigsaw::IterativeSinSaw>::new(),
        BandLimitedSignal::<ReferenceSaw>::new(),
        None,
        |bucket| bucket.data() == 0,
    )
    .unwrap()
}

#[test]
#[ignore]
/// Compare the saw with multiply-by-inverse to the reference saw
fn invmul_vs_reference_saw() {
    shared::logger::init_logger();
    plot_initial_error(
        BandLimitedSignal::<jigsaw::InvMulSaw>::new(),
        BandLimitedSignal::<ReferenceSaw>::new(),
        None,
        |bucket| bucket.data() == 0,
    )
    .unwrap()
}

#[test]
#[ignore]
/// Compare the saw using a smart FFT-like harmonics computation method
/// to the reference saw
fn smartharms_vs_reference_saw() {
    shared::logger::init_logger();
    plot_initial_error(
        BandLimitedSignal::<jigsaw::SmartHarmonicsSaw>::new(),
        BandLimitedSignal::<ReferenceSaw>::new(),
        None,
        |bucket| bucket.data() == 0,
    )
    .unwrap()
}

#[test]
#[ignore]
/// Compare the saw using a fully iterative harmonics computation method to the
/// reference saw
fn fullit_vs_reference_saw() {
    shared::logger::init_logger();
    plot_initial_error(
        BandLimitedSignal::<jigsaw::FullyIterativeSaw>::new(),
        BandLimitedSignal::<ReferenceSaw>::new(),
        None,
        |bucket| bucket.data() == 0,
    )
    .unwrap()
}
