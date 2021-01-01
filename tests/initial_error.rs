//! This test studies the error of band-limited signals at their initial phase,
//! both with respect to each other and to band-unlimited signals.

mod error;
mod logger;
mod map;
pub mod parameters;
pub mod signal;

use crate::{
    error::measure_error,
    logger::init_logger,
    map::OscillatorMap,
    parameters::{
        log2_relative_rate_range, NUM_PHASE_BUCKETS, NUM_RELATIVE_FREQ_BUCKETS, PHASE_RANGE,
    },
    signal::{BandLimitedSignal, Signal, UnlimitedSignal},
};
use core::sync::atomic::AtomicU8;
use jigsaw::{AudioPhaseMod, AudioSample, ReferenceSaw};
use log::{info, trace};
use plotters::prelude::*;

/// Compare two implementations of a given signal
fn plot_initial_error(
    signal: impl Signal + Send + Sync,
    reference: impl Signal + Send + Sync,
    plot_filename: &str,
) -> Result<(), Box<dyn std::error::Error>> {
    // Prepare to plot the error data
    use AudioPhaseMod::consts::PI;
    const X_MARGIN: u32 = 50;
    const Y_MARGIN: u32 = 60;
    const LABEL_SIZE: u32 = 20;
    const NUM_X_LABELS: usize = 12;
    const NUM_Y_LABELS: usize = 18;
    let relative_freq_buckets = NUM_RELATIVE_FREQ_BUCKETS as u32;
    let phase_buckets = NUM_PHASE_BUCKETS as u32;
    let phase_range_in_pi = (PHASE_RANGE.start / PI)..(PHASE_RANGE.end / PI);
    let root = BitMapBackend::new(
        plot_filename,
        (relative_freq_buckets + Y_MARGIN, phase_buckets + X_MARGIN),
    )
    .into_drawing_area();
    root.fill(&WHITE)?;
    let mut chart = ChartBuilder::on(&root)
        .margin(0)
        .x_label_area_size(X_MARGIN)
        .y_label_area_size(Y_MARGIN)
        .build_cartesian_2d(log2_relative_rate_range(), phase_range_in_pi)?;
    chart
        .configure_mesh()
        .disable_mesh()
        .label_style(("sans-serif", LABEL_SIZE))
        .x_desc("Sampling rate / oscillator frequency")
        .x_label_formatter(&|&log2_x| format!("{:.0}", (2.0f32).powf(log2_x)))
        .x_labels(NUM_X_LABELS)
        .y_desc("Phase in units of pi")
        .y_labels(NUM_Y_LABELS)
        .draw()?;
    let plotting_area = chart.plotting_area();
    assert_eq!(
        plotting_area.dim_in_pixel(),
        (relative_freq_buckets, phase_buckets)
    );
    let (base_x, base_y) = plotting_area.get_base_pixel();

    // Map the error landscape
    let error_map = OscillatorMap::measure(
        |sampling_rate, oscillator_freq, initial_phase| {
            measure_error(
                &signal,
                &reference,
                sampling_rate,
                oscillator_freq,
                initial_phase,
            )
        },
        AtomicU8::fetch_max,
    );

    // Draw the error map
    trace!("Error(relfreq, phase) table is:");
    trace!("relrate,phase,error");
    let mut max_error = 0;
    for bucket in error_map.iter() {
        // Check out the current error map bucket
        // FIXME: For this to become a real test, we should also inspect the error
        let error = bucket.data();
        max_error = max_error.max(error);
        let (relative_rate, phase) = bucket.center();
        trace!("{},{},{}", relative_rate, phase, error);

        // Plot the current bucket on the error map
        let scaled_error = ((error as f32 / AudioSample::MANTISSA_DIGITS as f32) * 255.0) as u8;
        let (x_idx, y_idx) = bucket.index();
        root.draw_pixel(
            (x_idx as i32 + base_x, y_idx as i32 + base_y),
            &RGBColor(scaled_error, 0, 0),
        )?;
    }
    info!("Maximum error is {} bits", max_error);
    Ok(())
}

#[test]
#[ignore]
/// Compare the reference saw to a band-unlimited saw
fn reference_vs_unlimited_saw() {
    init_logger();
    plot_initial_error(
        BandLimitedSignal::<ReferenceSaw>::new(),
        UnlimitedSignal::new(jigsaw::unlimited_saw),
        "initial_reference_vs_unlimited_saw.png",
    )
    .unwrap()
}

#[test]
#[ignore]
/// Compare the saw with single-precision sinus to the reference saw
fn f32sin_vs_reference_saw() {
    init_logger();
    plot_initial_error(
        BandLimitedSignal::<jigsaw::F32SinSaw>::new(),
        BandLimitedSignal::<ReferenceSaw>::new(),
        "initial_f32sin_vs_reference_saw.png",
    )
    .unwrap()
}

#[test]
#[ignore]
/// Compare the saw with iterative sinus to the reference saw
fn itersin_vs_reference_saw() {
    init_logger();
    plot_initial_error(
        BandLimitedSignal::<jigsaw::IterativeSinSaw>::new(),
        BandLimitedSignal::<ReferenceSaw>::new(),
        "initial_itersin_vs_reference_saw.png",
    )
    .unwrap()
}

#[test]
#[ignore]
/// Compare the saw with multiply-by-inverse to the reference saw
fn invmul_vs_reference_saw() {
    init_logger();
    plot_initial_error(
        BandLimitedSignal::<jigsaw::InvMulSaw>::new(),
        BandLimitedSignal::<ReferenceSaw>::new(),
        "initial_invmul_vs_reference_saw.png",
    )
    .unwrap()
}

#[test]
#[ignore]
/// Compare the saw using a smart FFT-like harmonics computation method
/// to the reference saw
fn smartharms_vs_reference_saw() {
    init_logger();
    plot_initial_error(
        BandLimitedSignal::<jigsaw::SmartHarmonicsSaw>::new(),
        BandLimitedSignal::<ReferenceSaw>::new(),
        "initial_smartharms_vs_reference_saw.png",
    )
    .unwrap()
}

#[test]
#[ignore]
/// Compare the saw using a fully iterative harmonics computation method to the
/// reference saw
fn fullit_vs_reference_saw() {
    init_logger();
    plot_initial_error(
        BandLimitedSignal::<jigsaw::FullyIterativeSaw>::new(),
        BandLimitedSignal::<ReferenceSaw>::new(),
        "initial_fullit_vs_reference_saw.png",
    )
    .unwrap()
}
