//! This test measures the error of band-limited signals with respect to their
//! non-band-limited equivalent through Monte-Carlo integration.
//!
//! Said error is a mixture of essential error from the band limiting and
//! accidental error from approximations in the computation. As the number of
//! harmonics in the band-limited signal goes up, both of these contributions go
//! up, but they do so differently:
//!
//! - In continuous portions of the input function, the error goes down as the
//!   number of harmonics increases. But for a typical function with Fourier
//!   coefficients of order 1/n, the precision improvement by adding the n-th
//!   harmonic is of order 1/n, so there are diminishing returns.
//!     * If we look at things from a geometric point of view, then a Fourier
//!       series that's truncated to n terms is of order of the sum of harmonic
//!       numbers sum(i in 1..=n, N(i)) = (n + 1) * H(n) - n. If we double
//!       the amount of terms, we get instead sum(i in 1..=2*n, N(i)), which is
//!       (2*n + 1) * H(2*n) - 2*n. So the difference in magnitude is of order
//!       (2*n + 1) * H(2*n) - 2*n - (n + 1) * H(n) + n
//!         = (2*n + 1) * H(2*n) - n - (n + 1) * H(n)
//!     * At large n, H(n) converges to gamma + ln(n), so we get approximately
//!       (2*n + 1) * (gamma + ln(2*n)) - n - (n + 1) * (gamma + ln(n))
//!         = (2*n + 1) * (gamma + ln(2) + ln(n)) - n - (n + 1) * (gamma + ln(n))
//!       The non-constant terms dominate the other ones, so we can simplify to
//!       (2*n + 1) * ln(n) - n - (n + 1) * ln(n)
//!         = n * ln(n) - n
//!         = n * (ln(n) - 1)
//!     * In other words, doubling the number of Fourier series terms should
//!       only give us a roughly linear improvement to the approximation error.
//!       Of course, that's for functions with infinite Fourier series like
//!       sawtooths, things like additive signal synthesis of real-world signals
//!       should eventually yiels a perfectly accurate clone of the original,
//!       and since we're computing with floating-point numbers, we'll sooner or
//!       later reach the point where adding extra harmonics doesn't improve the
//!       integration error at all.
//! - Around discontinuities, the Gibbs phenomenon error increases as the amount
//!   of harmonics increases, and converges to an error of a bit less than 9% on
//!   both sides of the discontinuity. But Gibbs peaks also become more and more
//!   narrow, which means that the odds of observing a significant trace of the
//!   Gibbs peaks on any sample of the signal goes down. Concretely, regions of
//!   the signal that are subjected to the Gibbs phenomenon should shrink and
//!   eventually become unobservable given enough harmonics.
//! - Computation error is a combination of harmonics computation error and
//!   accumulation error. It is localized in areas where either of these
//!   contributions is maximal.
//!
//! Assuming that the libm sine function is perfectly precise (within 1/2 ulp of
//! the true mathematical sine), the computation error of the reference sine
//! that computes and accumulates sines using double precision can be estimated:
//!
//! - Harmonics count are not allowed to overflow an f32 mantissa, so the
//!   floating-point error of casting the harmonic counter to f64 is zero.
//! - Phase conversion from f32 to f64 is also lossless.
//! - Multiplying the harmonics by the phase yields an error of 1/2 f64 ulp.
//! - Sine is between 0 and 1, so if the sine is computed within 1/2 ulp of the
//!   true mathematical result, the error on sine(n*phi +/- 1/2 ulp) is
//!   less than 1 f64 epsilon (1/2 eps from the sine computation, and 1/2 ulp
//!   from the input that propagate into the sine whose derivative is at most 1
//!   and whose result is at most 1, yielding a propagated error < 1/2 eps)
//! - Dividing the sine by the harmonic counter yields an error of 1/2 f64 ulp,
//!   on a result whose order of magnitude is 1/n where n is the harmnonic
//!   counter, and it scales the sine error by 1/n. So we get (1/2n) eps error
//!   from the division and at worst 1/n eps error from the sine, for a total
//!   error of at worst (3/2n) eps on a result whose order of magnitude is 1/n.
//!   In other words, the error on sin(nx) / n is at worst 3/2 ulp.
//! - Now, let's study the accumulation error, first assuming that we sum from
//!   the first harmonic to the last harmonic:
//!     * Intuitively, this error should be similar to that of computing an
//!       harmonic number N(i) = sum(j in 1..=i, 1/j), which itself should be
//!       proportional to said harmonic number because we're doing float
//!       computations and the precision of basic float operations is relative
//!       to the magnitude of their result due to round-off. So the final
//!       expression should feature harmonic number terms, and we'll try to make
//!       those emerge as needed. For that, we'll need the recurence relation
//!       N(i+1) = (N(i) + 1/i).
//!     * Also, because the error of basic float operations is 1/2 ulp, we'll
//!       expect something that's proportional to 1/2, so we'll try to extract
//!       1/2 prefactors as needed.
//!     * On the first sum, we sum something of order 1 with an error of order
//!       (1/2 * 3) eps and something of order 1/2 with an error of order
//!       (1/2 * 3 * 1/2) eps.
//!       If the sum is reasonably accurate, the result is of order N(2),
//!       and since sums have an error of 1/2 ulp that gives us a summation
//!       error of order [1/2 * N(2)] eps.
//!       The total error, when combining term errors and sum error, is thus:
//!       {1/2 * [3 + 1/2 * 3 + N(2)]} eps
//!         = {1/2 * [3 * N(2) + N(2)]} eps
//!         = [1/2 * 4 * N(2)] eps
//!       The relative error is of order of that divided by N(2), i.e.
//!       [1/2 * 4 * N(2) / N(2)] ulp
//!         = (1/2 * 4) ulp
//!         = 2 ulp
//!     * On the second sum, we sum an accumulator of order N(2) with an error
//!       error of [1/2 * 4 * N(2)] eps with a new term of order 1/3 with an
//!       error of order [1/2 * 3 * 1/3] eps.
//!       If the sum is reasonably accurate, the result is of order
//!       N(3), and since sums have an error of 1/2 ulp that gives us a
//!       summation error of order [1/2 * N(3)] eps.
//!       The total error, when combining term errors and sum error, is thus:
//!       {1/2 * [4 * N(2) + 3 * 1/3 + N(3)]} eps
//!         = {1/2 * [N(2) + 3 * (N(2) + 1/3) + N(3)]} eps
//!         = {1/2 * [N(2) + 3 * N(3) + N(3)]} eps
//!         = {1/2 * [N(2) + 4 * N(3)]} eps
//!       The relative error is of order of that divided by N(3), i.e.
//!       {1/2 * [N(2) + 4 * N(3)] / N(3)} ulp
//!         = {1/2 * [N(2) / N(3) + 4]} ulp
//!         = [2 + 1/2 * N(2) / N(3)] ulp
//!     * On the third sum, we sum an accumulator of order N(3) with an error of
//!       {1/2 * [N(2) + 4 * N(3)]} eps with a new term of order 1/4 with an
//!       error of order [1/2 * 3 * 1/4] eps.
//!       If the sum is reasonably accurate, the result is of order
//!       N(4), and since sums have an error of 1/2 ulp that gives us a
//!       summation error of order [1/2 * N(4)] eps.
//!       The total error, when combining term errors and sum error, is thus:
//!       {1/2 * [N(2) + 4 * N(3) + 3 * 1/4 + N(4)]} eps
//!         = {1/2 * [N(2) + N(3) + 3 * (N(3) + 1/4) + N(4)]} eps
//!         = {1/2 * [N(2) + N(3) + 3 * N(4) + N(4)]} eps
//!         = {1/2 * [N(2) + N(3) + 4 * N(4)]} eps
//!       The relative error is of order of that divided by N(4), i.e.
//!       {1/2 * [N(2) + N(3) + 4 * N(4)] / N(4)} ulp
//!         = {2 + 1/2 * [N(2) + N(3)] / N(4)} ulp
//!     * On the fourth sum, we sum an accumulator of order N(4) with an error
//!       of {1/2 * [N(2) + N(3) + 4 * N(4)]} eps with a new term of order 1/5
//!       with an error of order [1/2 * 3 * 1/5] eps.
//!       If the sum is reasonably accurate, the result is of order N(5), and
//!       since sums have an error of 1/2 ulp that gives us a summation error of
//!       order [1/2 * N(5)] eps.
//!       The total error, when combining term errors and sum error, is thus:
//!       {1/2 * [N(2) + N(3) + 4 * N(4) + 3 * 1/5 + N(5)]} eps
//!         = {1/2 * [N(2) + N(3) + N(4) + 3 * (N(4) + 1/5) + N(5)]} eps
//!         = {1/2 * [N(2) + N(3) + N(4) + 3 * N(5) + N(5)]} eps
//!         = {1/2 * [N(2) + N(3) + N(4) + 4 * N(5)]} eps
//!       The relative error is of order of that divided by N(5), i.e.
//!       {1/2 * [N(2) + N(3) + N(4) + 4 * N(5)] / N(5)} ulp
//!         = {2 + 1/2 * [N(2) + N(3) + N(4)] / N(5)} ulp
//!     * At this point, it's pretty clear that in the whole range where the
//!       sum is reasonably accurate, the absolute error of the i-th sum is
//!       very likely to be
//!       {1/2 * [sum(j in 2..=i, N(j)) + 3 * N(i)]} eps
//!         = [1/2 * sum(j in 2..=i, N(j)) + 3/2 * N(i)] eps
//!         = [-1/2 + 1/2 * sum(j in 1..=i, N(j)) + 3/2 * N(i)] eps
//!     * People who spent more time than I'm willing to studying harmonic
//!       series tell us that sum(j in 1..=i, N(j)) = (i + 1) * N(i) - i
//!       Therefore the absoute error is actually
//!       {-1/2 + 1/2 * [(i + 1) * N(i) - i] + 3/2 * N(i)} eps
//!         = {-1/2 + 1/2 * [(i + 4) * N(i) - i]} eps
//!         = {1/2 * [-1 + (i + 4) * N(i) - i]} eps
//!     * At large i, the constant terms become negligible, and thus the error
//!       becomes of the order of
//!       {1/2 * [i * N(i) - i]} eps
//!         = {1/2 * i * (N(i) - 1)} eps
//!       And by the same trick, we get rid of that -1 at very large i, which
//!       gives us an error of around
//!       {1/2 * i * N(i)} eps
//!       And thus a relative error of around i/2 ulp.
//!     * This means that if we use a double-precision accumulator, we'll only
//!       see an observable error when casting to single-precision if we
//!       accumulate more terms than 2^(s+1) where s is the number of bits in
//!       the mantissa of a single-precision float. OTOH, doing the whole
//!       accumulation in single precision would clearly be a very bad idea, at
//!       least for a large number of harmonics.
//!     * Anyway, what this tells us is that for the reference saw, which does
//!       everything in double precision and casts to single precision at the
//!       end, all observable error should be attributable to missing Fourier
//!       terms and the Gibbs phenomenon, which we analyzed above.
//!
//! Now that the maths have been spelled out, let's check that our reference
//! band-limited saw implementation does honor them.
//!
//! To do that, we'll estimate the shape of the error function representing the
//! difference between the band-limited and unlimited parameters through Monte
//! Carlo integration. More concretely...
//!
//! - We loop over (sampling rate, freq) in the right range (audio frequencies,
//!   sampling rates between 44100 Hz and 2^f32::MANTISSA_BITS).
//!     * We generate uniform random phases in the 0..2*pi range.
//!         - For each phase, we estimate the integration error, as measured by
//!           the figure of merit of the number of wrong bits
//!           ((signal1 - signal2).log2().ceil()).
//!         - We do that by computing the error on successive samples of the
//!           band-limited signal and their band-limited equivalent, taking the
//!           max of those errors until said max saturates for N successive
//!           samples (this reduces noise and amortizes the overhead of
//!           constructing a signal a bit, at the cost of some loss of
//!           spatial resolution in the phase domain.
//!         - With that, we keep updated a list of "error domains" where the
//!           error is bounded by 2^0 bits, 2^1 bits, 2^2 bits...
//!         - We stop integrating when the number of error domains converges (it
//!           has not increased for M1 successive samples).
//!         - We can then refine domain boundaries by sampling more data points
//!           in the regions between domains. If we see more domains by doing
//!           so, it means that M1 is not big enough, so we abort with an error.
//!     * We repeat the process for various (sampling rate, freq) configurations
//!       and we store the results by relative frequency as that should be the
//!       right parameter for what we are studying.
//!     * Our stopping criteria is to look at rectangular "error domains" in the
//!       transverse direction, and again stop when their number stops
//!       increasing for a number of iterations M2.
//!     * Then we can refine error domains by looking at relative frequencies
//!       between the domains, again aborting with an error if we find new
//!       domains while doing so, as it means that M2 is too small.
//! - At the end, we output a list of error domain quads, which can be used to
//!   draw a map of the error vs relative frequency and phase. We could actually
//!   draw that map ourselves using something like plotters.

use jigsaw::{
    unlimited_saw, AudioFrequency, AudioPhase, AudioPhaseMod, AudioSample, Oscillator,
    OscillatorPhase, ReferenceSaw, SamplingRateHz,
};
use log::{debug, info, trace};
use rand::Rng;
use std::{
    cmp::Ordering,
    collections::VecDeque,
    ops::{Range, RangeInclusive},
};

// Sampled region of the oscillator configuration space
const SAMPLING_RATE_RANGE: RangeInclusive<SamplingRateHz> =
    44100..=(1 << AudioFrequency::MANTISSA_DIGITS);
const OSCILLATOR_FREQ_RANGE: RangeInclusive<AudioFrequency> = 20.0..=20000.0;
const PHASE_RANGE: Range<AudioPhase> = 0.0..AudioPhaseMod::consts::TAU;

/// Phase range over which error should be integrated
///
/// An ideal Monte Carlo sampler only takes one sample, but there are several
/// benefits to taking more of them per considered oscillator phases:
///
/// - It amortizes the overhead of constructing an oscillator, and allows
///   runtime performance optimizations like vectorization to kick in.
/// - It produces a smoothed out error domain landscape that will converge
///   faster and use less RAM, but may still provide enough information.
///
const ERROR_INTEGRATION_RANGE: AudioPhase = AudioPhaseMod::consts::TAU / 1000.0;

/// Number of bits of error from an unlimited signal to its limited equivalent
type ErrorBits = u8;

/// Error integration data point
type ErrorDataPoint = (AudioPhase, ErrorBits);

/// Measure how the reference saw differs from the band-unlimited saw
///
/// Returns the central phase around which the error was measured, and the
/// number of bits that differ as a coarse figure of merit of the error.
///
fn measure_error(
    sampling_rate: SamplingRateHz,
    oscillator_freq: AudioFrequency,
    initial_phase: AudioPhase,
) -> ErrorDataPoint {
    let oscillator = ReferenceSaw::new(sampling_rate, oscillator_freq, initial_phase);
    let osc_phase = OscillatorPhase::new(sampling_rate, oscillator_freq, initial_phase);
    debug!(
        "Probing error at phase {} (= {}pi)",
        initial_phase,
        initial_phase / AudioPhaseMod::consts::PI
    );
    let mut integrated_error_bits = 0;
    let mut final_phase = 0.0;
    trace!("Phase,BandLimited,TrueSaw,Difference,ErrorBits,IntegratedErrorBits");
    for (phase, limited) in osc_phase.zip(oscillator) {
        use AudioPhaseMod::consts::TAU;
        let unlimited = unlimited_saw(phase);
        let difference = limited - unlimited;
        let correct_bits = (-(difference.abs().log2()))
            .floor()
            .max(0.0)
            .min(u8::MAX as AudioSample) as u8;
        let error_bits = (AudioSample::MANTISSA_DIGITS as u8).saturating_sub(correct_bits);
        integrated_error_bits = integrated_error_bits.max(error_bits);
        trace!(
            "{},{},{},{},{},{}",
            phase,
            limited,
            unlimited,
            difference,
            error_bits,
            integrated_error_bits
        );
        if (phase - initial_phase).rem_euclid(TAU) >= ERROR_INTEGRATION_RANGE {
            final_phase = phase.rem_euclid(TAU);
            break;
        }
    }
    let average_phase = (initial_phase + final_phase) / 2.0;
    debug!(
        "Went up to phase {} (= {}pi), found error of {} bits",
        final_phase,
        final_phase / AudioPhaseMod::consts::PI,
        integrated_error_bits
    );
    (average_phase, integrated_error_bits)
}

/// From these error measurements, we integrate an increasingly accurate model
/// of the error landscape of the function.
///
/// To enable reasonably fast convergence, this model is limited to a spatial
/// resolution of ERROR_INTEGRATION_RANGE. Concretely, if an error data point
/// has a neighbour with that range which has a higher error, it is considered
/// to have that error as well.
///
/// This leads to the iterative formation and stabilization of "error domains"
/// which are pessimistically considered to have constant integration error.
///
// TODO: Spell out the full list of invariants and test them in debug mode
struct ErrorLandscapeIntegrator {
    /// List of error domains, ordered by min/max phase
    domain_integrators: Vec<ErrorDomainIntegrator>,
}
//
struct ErrorDomainIntegrator {
    /// Upper bound on the number of bits of error within that domain
    error: ErrorBits,

    /// Ordered list of phases that were associated to that error domain
    ///
    /// There should always be at least one phase data point within an error
    /// domain. If a domain falls to zero phases, it should be destroyed.
    ///
    phases: VecDeque<ErrorDataPoint>,
}
//
struct ErrorDomain {
    /// Phase at which the domain starts
    start: AudioPhase,

    /// Number of bits of error within that domain
    error: ErrorBits,

    /// Phase at which the domain ends
    end: AudioPhase,
}
//
impl ErrorLandscapeIntegrator {
    /// Prepare to integrate the error landscape
    fn new() -> Self {
        Self {
            domain_integrators: Vec::new(),
        }
    }

    /// Query how many domains are present in the current error landscape
    fn num_domains(&self) -> usize {
        self.domain_integrators.len()
    }

    /// Query a domain within the current error landscape
    fn domain(&self, index: usize) -> ErrorDomain {
        self.domain_integrators[index].domain()
    }

    /// Query the current error domain landscape
    ///
    /// This iterator returns a list of domain, where for each domain, we get
    /// the phase at which the domain starts, the number of error bits within
    /// that error domain, and the phase at which the domain ends.
    ///
    fn domains(&self) -> impl Iterator<Item = ErrorDomain> + '_ {
        self.domain_integrators
            .iter()
            .map(ErrorDomainIntegrator::domain)
    }

    /// Update the error landscape with a new measurement
    ///
    /// Tell whether error domains were created or destroyed in this process, or
    /// the boundary of existing error domains was simply moved around.
    fn integrate(&mut self, data_point: ErrorDataPoint) -> bool {
        // Separate the components of the error measurement
        let (average_phase, error) = data_point;

        // Handle a measurement that falls within an existing error domain
        let integrate_into = |integrators: &mut Vec<ErrorDomainIntegrator>, domain_idx| {
            let domain_integrator: &mut ErrorDomainIntegrator = &mut integrators[domain_idx];
            let domain_error = domain_integrator.error;
            match domain_error.cmp(&error) {
                // Same error: Just insert the data point there
                Ordering::Equal => {
                    domain_integrator.insert_data(data_point);
                    false
                }

                // FIXME: Find out how to avoid high-error domain gradually
                //        absorbing all data through percolation. I think
                //        keeping track of the actual error within each data
                //        point of a domain is key, which is why I now record
                //        this info, but algorithm must also change to use that.
                //
                // Domain has more error than input data point and would tend to
                // absorb said input, but will only do so if it contains a close
                // enough phase. Otherwise, the domain will be split and a new
                // phase domain will be inserted in the middle.
                Ordering::Greater => {
                    /* trace!("Candidate has more error bits, will it absorb us?");
                    match search_phase(candidate, average_phase, true) {
                        Ok(_approx_pos) => {
                            trace!("Yes, it will absorb us");
                            insert_phase(candidate)
                        }
                        Err(_approx_splitting_point) => {
                            trace!("No, we split it and insert ourselves in the middle");
                            // FIXME: Handle case where splitting point is zero
                            let right_candidate_phases =
                                candidate.phases.split_off(splitting_point);
                            error_domains.insert(
                                match_idx + 1,
                                ErrorDomain {
                                    phases: std::iter::once(average_phase).collect(),
                                    error_bits: error_bits,
                                },
                            );
                            error_domains.insert(
                                match_idx + 2,
                                ErrorDomain {
                                    phases: right_candidate_phases,
                                    error_bits: error_bits,
                                },
                            );
                        }
                    } */
                    todo!()
                }

                // FIXME: Find out how to avoid high-error domain gradually
                //        absorbing all data through percolation. I think
                //        keeping track of the actual error within each data
                //        point of a domain is key, which is why I now record
                //        this info, but algorithm must also change to use that.
                //
                // Domain has less error than input data point and all of its
                // data points that lie close to the input phase will be
                // absorbed by a new domain created for the input data point.
                Ordering::Less => {
                    /*
                    trace!("Candidate has less error bits, create new domain and absorb its closest phases");
                    let left_splitting_point = search_phase(
                        candidate,
                        average_phase - ERROR_INTEGRATION_RANGE,
                        false,
                    )
                    .unwrap_or_else(|pos| pos);
                    // TODO: Split candidate in three parts:
                    // - Phases that are far below us (left_candidate_phases)
                    // - Phases that are close to us (close_candidate_phases)
                    // - Phases that are far above us (right_candidate_phases)
                    // If left_candidate_phases is not empty, it stays where the candidate
                    // was. Then we insert a domain with our phase, close_candidate_phases,
                    // and our integrated error bits. Finally, we insert a domain with
                    // left_candidate_phases and candidate_error_bits.*/
                    todo!()
                }
            }
        };

        // Find if input phase falls inside a domain or between two domains
        match self.domain_integrators.binary_search_by(|probe| {
            let domain = probe.domain();
            if domain.start > average_phase {
                Ordering::Greater
            } else if domain.end < average_phase {
                Ordering::Less
            } else {
                Ordering::Equal
            }
        }) {
            // Input phase lies inside of an existing error domain
            Ok(idx) => {
                trace!("Integrating into existing error domain {}", idx);
                integrate_into(&mut self.domain_integrators, idx)
            }

            // The phase lies outside existing error domains
            Err(idx) => {
                // Check if neighboring domains are within integration range or
                // have the same error value, which means that the active
                // measurement will interact with them.
                let left_border = if idx > 0 {
                    let left_domain = self.domain(idx - 1);
                    average_phase - ERROR_INTEGRATION_RANGE <= left_domain.end
                        || left_domain.error == error
                } else {
                    false
                };
                let right_border = if idx < self.num_domains() {
                    let right_domain = self.domain(idx);
                    average_phase + ERROR_INTEGRATION_RANGE >= right_domain.end
                        || right_domain.error == error
                } else {
                    false
                };

                // Deduce what action needs to be taken
                match (left_border, right_border) {
                    // Input measurement does not match any neighboring ones, so
                    // we can create a new error domain just for it
                    (false, false) => {
                        self.domain_integrators.insert(
                            idx,
                            ErrorDomainIntegrator {
                                error,
                                phases: std::iter::once(data_point).collect(),
                            },
                        );
                        true
                    }

                    // Input measurement matches one neighboring domain and
                    // effectively belongs to that domain as in the case above
                    (true, false) => integrate_into(&mut self.domain_integrators, idx - 1),
                    (false, true) => integrate_into(&mut self.domain_integrators, idx),

                    // Phase belongs to two domains at idx-1 and idx, must arbitrate
                    (true, true) => {
                        // FIXME: Find out how to avoid high-error domain gradually
                        //        absorbing all data through percolation. I think
                        //        keeping track of the actual error within each data
                        //        point of a domain is key, which is why I now record
                        //        this info, but algorithm must also change to use that.
                        //
                        // TODO: By construction, this can only happen when the
                        //       phase lies inbetween two domains, which should
                        //       simplify matters. However, we can still "build
                        //       a bridge" between two domains, which means that
                        //       the action that we need to take is more
                        //       complicated than merging into the first domain,
                        //       then the second domain.
                        todo!()
                    }
                }
            }
        };
        todo!("Finish this")
    }
}
//
impl ErrorDomainIntegrator {
    /// Extract the current domain statistics for this domain
    fn domain(&self) -> ErrorDomain {
        ErrorDomain {
            start: self.phases.front().unwrap().0,
            error: self.error,
            end: self.phases.back().unwrap().0,
        }
    }

    /// Search function for finding a certain phase in an error
    // domain through binary search, with or without a tolerance
    fn search_phase(&mut self, phase: AudioPhase, with_tolerance: bool) -> Result<usize, usize> {
        self.phases
            .make_contiguous()
            .binary_search_by(|&(probed_phase, _probed_error)| {
                let tolerance = (with_tolerance as u8 as AudioPhase) * ERROR_INTEGRATION_RANGE;
                if probed_phase > phase + tolerance {
                    Ordering::Greater
                } else if probed_phase < phase - tolerance {
                    Ordering::Less
                } else {
                    Ordering::Equal
                }
            })
    }

    // Insert a data point at the right position within an error domain
    fn insert_data(&mut self, data_point: ErrorDataPoint) {
        match self.search_phase(data_point.0, false) {
            Err(pos) => {
                trace!("Phase wasn't previously probed, record it");
                self.phases.insert(pos, data_point);
            }
            Ok(pos) => {
                trace!("Phase was already probed before");
                debug_assert_eq!(self.phases[pos], data_point);
            }
        }
    }
}

#[test]
#[ignore]
fn reference_saw() {
    env_logger::init();
    let mut rng = rand::thread_rng();
    loop {
        // Pick a relative oscillator frequency
        let sampling_rate = rng.gen_range(SAMPLING_RATE_RANGE);
        let oscillator_freq = rng.gen_range(OSCILLATOR_FREQ_RANGE);
        info!(
            "Sampling a {} Hz saw at {} Hz (relative freq = {})",
            oscillator_freq,
            sampling_rate,
            oscillator_freq / sampling_rate as AudioFrequency
        );

        // Map the error landscape at that relative frequency
        let mut error_landscape = ErrorLandscapeIntegrator::new();
        loop {
            // Check error at random phases
            let initial_phase = rng.gen_range(PHASE_RANGE);
            let error_measurement = measure_error(sampling_rate, oscillator_freq, initial_phase);

            // Use these measurements to study the error landscape via a form
            // of Monte Carlo integration
            if error_landscape.integrate(error_measurement) {
                debug!(
                    "Error landscape changed, there are now {} error domains",
                    error_landscape.domains().count()
                );
            }

            // TODO: Stop when error domain list converges. But currently, it
            //       doesn't converge. So I must first figure out why.
            if error_landscape.domains().count() > 10000 {
                debug!("Index,StartPhase,EndPhase,ErrorBits");
                for (idx, domain) in error_landscape.domains().enumerate() {
                    debug!("{},{},{},{}", idx, domain.start, domain.end, domain.error,);
                }
                todo!()
            }
        }
        // TODO: Do something useful
        todo!()
    }
    // TODO: Do something useful
    todo!()
}
