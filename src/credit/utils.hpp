#ifndef UTILS_H
#define UTILS_H

#include <cmath>
#include <concepts>
#include <cstdint>
#include <stdexcept>

#include <boost/math/tools/toms748_solve.hpp>

inline constexpr double TOL = 1e-3;
inline constexpr int N_ITER = 50;

/* Description
 * -----------
 * Root solver via a guaranteed-convergence bracketing method.
 *
 * Why a bracketing method (and not the previous secant iteration):
 * - The secant update divides by f(x2) - f(x1). Our objective functions contain
 *   Black-Scholes prices, which are numerically flat deep in/out of the money, so
 *   two nearby evaluations can be equal to double precision. That makes the secant
 *   step divide by zero (or take an enormous step from a near-zero denominator).
 * - A bracketing method instead maintains an interval [a, b] with f(a) and f(b) of
 *   opposite sign. A root is then guaranteed to lie inside, every iteration shrinks
 *   the interval, and no step can diverge or divide by zero: worst case it degrades
 *   to bisection, which always converges.
 * - Brent's method is the classic algorithm of this family (bisection + secant +
 *   inverse quadratic interpolation). Boost.Math ships TOMS 748, its direct
 *   successor: the same bracketing guarantee, but with cubic instead of quadratic
 *   interpolation, so it typically needs fewer function evaluations. We use that
 *   rather than hand-rolling Brent.
 *
 * The bracket starts at [x0/2, 2*x0] (the neighborhood the old secant solver used)
 * and expands geometrically until a sign change is found. Halving/doubling keeps the
 * bracket strictly positive, which matters because callers probe option prices where
 * the underlying value must stay > 0.
 *
 * Params
 * ------
 *  Func &&f:           Function whose root is sought
 *  T x0:               Initial guess; must be positive
 *  T tol = 1e-3:       Absolute width of the final bracket
 *  int n_iter = 50:    Maximum number of refinement iterations (and expansions)
 *
 * Returns
 * -------
 * T x:                Root solution (midpoint of the final bracket)
 *
 * Throws std::invalid_argument if x0 <= 0, and std::runtime_error if no sign
 * change can be bracketed, i.e. the function never crosses zero in reach.
 */
template <typename Func, std::floating_point T>
    requires std::invocable<Func&, T>
T bracket_root(Func&& f, T x0, T tol = TOL, int n_iter = N_ITER)
{
    if (x0 <= T(0))
        throw std::invalid_argument("bracket_root: initial guess must be positive");

    auto a = x0 / 2;
    auto b = x0 * 2;
    auto fa = f(a);
    auto fb = f(b);

    // Expand the bracket geometrically until f changes sign across it
    for (auto iter = 0; iter < n_iter && fa * fb > T(0); ++iter)
    {
        a /= 2;
        b *= 2;
        fa = f(a);
        fb = f(b);
    }

    // Either endpoint may already be a root
    if (fa == T(0)) return a;
    if (fb == T(0)) return b;

    if (fa * fb > T(0))
        throw std::runtime_error(
            "bracket_root: could not bracket a sign change; function may have no root");

    auto stop = [tol](T lo, T hi) { return std::abs(hi - lo) <= tol; };
    std::uintmax_t max_iter = static_cast<std::uintmax_t>(n_iter);

    auto [lo, hi] = boost::math::tools::toms748_solve(f, a, b, fa, fb, stop, max_iter);

    return (lo + hi) / 2;
}

#endif
