//! Foundation for manipulating audio signals

/// Integer type suitable for storing an audio sampling rate in Hz
//
// In an age where high-end audio cards go to 192 000 Hz, u16 is too little, but
// u32 seems just fine for all foreseeable future.
//
pub type SamplingRateHz = u32;
#[cfg(test)]
pub type NonZeroSamplingRate = core::num::NonZeroU32;

/// Floating-point type suitable for storing audio frequencies
//
// Single precision enables storing frequencies with ~10^-7 precision, this is
// better than the human ear's frequency discrimination ability.
//
pub type AudioFrequency = f32;

/// Floating-point type suitable for storing audio samples
//
// Single precision gives us a precision equivalent to 24-bit integers, which
// are currently the gold standard of integer audio sample storage, and way
// beyond human earing's amplitude discrimination abilities. This is fine.
//
pub type AudioSample = f32;
