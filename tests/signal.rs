//! Common abstraction layer across band-limited and band-unlimited signals

use core::marker::PhantomData;
use jigsaw::{AudioFrequency, AudioPhase, AudioSample, Oscillator, SamplingRateHz};

/// Common trait to homogeneously handle the signal of band-limited oscillators
/// and that of the underlying band-unlimited mathematical function
pub trait Signal {
    /// Measure the value of signal at a certain phase, for certain parameters
    fn measure(
        &self,
        sampling_rate: SamplingRateHz,
        oscillator_freq: AudioFrequency,
        phase: AudioPhase,
    ) -> AudioSample;
}

/// Signal implementation based on an Oscillator implementation
pub struct BandLimitedSignal<Osc: Oscillator>(PhantomData<Osc>);
//
impl<Osc: Oscillator> BandLimitedSignal<Osc> {
    /// Implement the Signal abstraction using an Oscillator
    pub fn new() -> Self {
        Self(PhantomData)
    }
}
//
impl<Osc: Oscillator> Signal for BandLimitedSignal<Osc> {
    /// Measure the signal of an oscillator configured for a certain operating
    /// frequency and sampling rate, at a certain phase.
    fn measure(
        &self,
        sampling_rate: SamplingRateHz,
        oscillator_freq: AudioFrequency,
        phase: AudioPhase,
    ) -> AudioSample {
        Osc::new(sampling_rate, oscillator_freq, phase)
            .next()
            .unwrap()
    }
}

/// Signal implementation based on a band-unlimited mathematical function
pub struct UnlimitedSignal<Sig: Fn(AudioPhase) -> AudioSample>(Sig);
//
impl<Sig: Fn(AudioPhase) -> AudioSample> UnlimitedSignal<Sig> {
    /// Implement the Signal abstraction using a mathematical function
    pub fn new(signal: Sig) -> Self {
        Self(signal)
    }
}
//
impl<Sig: Fn(AudioPhase) -> AudioSample> Signal for UnlimitedSignal<Sig> {
    //// Measure the mathematical function at a certain phase
    fn measure(
        &self,
        _sampling_rate: SamplingRateHz,
        _oscillator_freq: AudioFrequency,
        phase: AudioPhase,
    ) -> AudioSample {
        (self.0)(phase)
    }
}
