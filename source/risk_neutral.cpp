// Implementation of prototypes in risk_neutral.h
#include <cmath>
using std::log;
using std::pow;

#include <boost/math/distributions/normal.hpp>
using boost::math::normal;

#include <boost/math/statistics/univariate_statistics.hpp>
// using boost::math::statistics::mean;
using boost::math::statistics::variance;

#include <boost/math/tools/roots.hpp>
using boost::math::tools::bisect;

#include "risk_neutral.h"

#ifndef TOL
#define TOL 1e-6
#endif

void option_price(
    const double S0,
    const double K,
    const double r,
    const double sigma,
    const double t,
    double& price,
    double& Phi1,
    double& Phi2,
    const bool call,
    const double q
) {
    double d1, d2;

    // Define N(0,1) object
    normal N(0.0, 1.0);

    d1 = (log(S0/K) + (r-q+0.5*pow(sigma, 2.0)) * t) / (sigma * sqrt(t));
    d2 = d1 - sigma * sqrt(t);

    if (call) {
        Phi1 = cdf(N, d1);
        Phi2 = cdf(N, d2);

        price = S0 * exp(-q*t) * Phi1 - K * exp(-r*t) * Phi2;
    } else {        
        Phi1 = cdf(N, -d1);
        Phi2 = cdf(N, -d2);

        price = K * exp(-r*t) * Phi2 - S0 * exp(-q*t) * Phi1;
    }
}

std::vector<double> wang_transform(
    std::vector<double> P,
    const double sharpe_ratio,
    const bool inverse
) {
    // Define N(0,1) object
    normal N(0.0, 1.0);

    for(auto& p : P) 
        p = cdf(N, quantile(N, p) + (inverse?-1:1) * sharpe_ratio);

    return P;
}

double _implied_volatility (
    double& e0,
    const double& S0,
    const double& K,
    const double& r,
    const double& sigma,
    const double& t
) {
    double impl_equity, Phi1, Phi2;
    option_price(S0, K, r, sigma, t, impl_equity, Phi1, Phi2);
    return (e0 - impl_equity);
}

double get_asset_volatility(
    const std::vector<double> E,
    const double sigma_e,
    const double L,
    const double r,
    const double t,
    unsigned int n_iter
) {
    // Check the trivial case: if no leverage, asset and equity volatility equate
    if (L < TOL)
        return sigma_e;

    // Obtain some properties of the equity price vector
    const unsigned int N = E.size();

    // Initialize guesses: A = E and sigma_a is the stdev of log returns
    std::vector<double> A(E);
    std::vector<double> log_returns(N-1);

    double sigma_a = 0.5;
    double prev_sigma_a;

    // Iterative steps
    for (unsigned int i = 0; i < n_iter; ++i) {
        // Save down a previous result for sigma_a
        prev_sigma_a = sigma_a; 

        // update sigma_a with current asset value vector
        for (unsigned int i = 1; i < N; ++i) 
            log_returns[i] = log(A[i] / A[i-1]);

        sigma_a = sqrt(variance(log_returns));

        // Using this guess of sigma_a, calculate the implied asset values in vector
        for (unsigned int j = 0; j < N; ++j) {
            const double& cur_equity = E[j];

            std::pair<double, double> result = bisect(
                // Solve for sigma_a to build out implied asset value
                [&cur_equity, &L, &r, &sigma_a, &t](double _a) { return _implied_volatility(_a, cur_equity, L, r, sigma_a, t); },
                0.0,
                1.0,
                [](double l, double r) { return abs(l-r) < TOL; }
            );

            // Save down new asset value
            A[j] = (result.first + result.second) / 2;
        }

        // Early stop if difference between previous and current sigma_iterations is below tolerance
        if (abs(sigma_a - prev_sigma_a) < TOL)
            break;
    }

    return sigma_a;
}

double get_default_probability(
    const double a0,
    const double mu_a,
    const double sigma_a,
    const double L,
    const double t
) {
    normal N(0.0, 1.0);

    double distance_to_default = (log(a0 / L) + (mu_a + 0.5 * pow(sigma_a, 2)) * t) / (sigma_a * sqrt(t));
    return cdf(N, -distance_to_default);
}
