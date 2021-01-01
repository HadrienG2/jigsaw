//! Common abstraction layer across band-limited and band-unlimited signals

use core::marker::PhantomData;
use jigsaw::{
    AudioFrequency, AudioPhase, AudioSample, Oscillator, OscillatorPhase, SamplingRateHz,
};

/// Common trait to homogeneously handle the signal of band-limited oscillators
/// and that of the underlying band-unlimited mathematical function
pub trait Signal {
    /// Iterator type which is produced when sampling this signal
    type Iter: Iterator<Item = AudioSample>;

    /// Start sampling a signal as an oscillator with certain parameters
    fn sample(
        &self,
        sampling_rate: SamplingRateHz,
        oscillator_freq: AudioFrequency,
        initial_phase: AudioPhase,
    ) -> Self::Iter;

    /// Measure the value of a signal at a certain phase
    fn measure(
        &self,
        sampling_rate: SamplingRateHz,
        oscillator_freq: AudioFrequency,
        phase: AudioPhase,
    ) -> AudioSample {
        self.sample(sampling_rate, oscillator_freq, phase)
            .next()
            .unwrap()
    }
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
    type Iter = Osc;

    fn sample(
        &self,
        sampling_rate: SamplingRateHz,
        oscillator_freq: AudioFrequency,
        initial_phase: AudioPhase,
    ) -> Self::Iter {
        Osc::new(sampling_rate, oscillator_freq, initial_phase)
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
impl<Sig: Fn(AudioPhase) -> AudioSample + Clone> Signal for UnlimitedSignal<Sig> {
    type Iter = std::iter::Map<OscillatorPhase, Sig>;

    fn sample(
        &self,
        sampling_rate: SamplingRateHz,
        oscillator_freq: AudioFrequency,
        initial_phase: AudioPhase,
    ) -> Self::Iter {
        OscillatorPhase::new(sampling_rate, oscillator_freq, initial_phase).map(self.0.clone())
    }
}
