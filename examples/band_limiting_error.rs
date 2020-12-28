//! This example measures the error of band-limited signals with respect to
//! their non-band-limited equivalent through Monte-Carlo integration.
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

fn main() {
    todo!()
}
