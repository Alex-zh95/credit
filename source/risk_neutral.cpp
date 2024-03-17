// Implementation of prototypes in risk_neutral.h
#include <cmath>
using std::log;
using std::pow;
using std::exp;
using std::sqrt;

#include <boost/math/distributions/normal.hpp>
using boost::math::normal;

#include "risk_neutral.hpp"

#ifndef TOL
#define TOL 1e-6
#endif

template <typename Func, typename T>
T secant_root(Func&& f, T& x0, T tol, int n_iter) {
    // Define a neighborhood around the initial value x0
    T x1 = x0/2;
    T x2 = x0*2;

    for (auto iter = 0; iter < n_iter; ++iter) {
        x0 = (x1 * f(x2) - x2 * f(x1)) / (f(x2) - f(x1));
        auto c = f(x0);

        if ((c < tol) && (c > -tol))
            break;

        if (c < -tol) {
            x1 = x2;
            x2 = x0;
        } else {
            x2 = x1;
            x1 = x0;
        }
    }

    return x0;
}

void vanilla_option_price(
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

double get_vanilla_asset_volatility(
    const std::vector<double> &E,
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
    std::vector<double> A = E;

    double sigma_a = 0.5;
    double prev_sigma_a;

    // Iterative steps
    for (unsigned int run_iter = 0; run_iter < n_iter; ++run_iter) {
        // Save down a previous result for sigma_a
        prev_sigma_a = sigma_a; 

        // update sigma_a with current asset value vector
        // From Itô Lemma, we have the following relationship: sigma_e * E[j] / A[j] = Phi1 * sigma_a
        double eq, Phi1, Phi2;
        vanilla_option_price(A[N-1], L, r, sigma_a, t, eq ,Phi1, Phi2);
        sigma_a = sigma_e * E[N-1] / A[N-1] / Phi1;

        // Early stop if difference between previous and current sigma_iterations is below tolerance
        if (abs(sigma_a - prev_sigma_a) < TOL)
            break;

        // Using this guess of sigma_a, calculate the implied asset values in vector
        for (unsigned int j = 0; j < N; ++j) {
            double cur_equity = E[j];

            auto fn = [&L, &r, &sigma_a, &t, &cur_equity](double _a) { 
                double impl_equity, Phi1, Phi2;
                vanilla_option_price(_a, L, r, sigma_a, t, impl_equity, Phi1, Phi2);
                return (cur_equity - impl_equity);
            };

            // New root is the asset value
            A[j] = secant_root(fn, cur_equity);
        }
    }

    return sigma_a;
}

double get_vanilla_default_probability(
    const double a0,
    const double mu_a,
    const double sigma_a,
    const double L,
    const double t
) {
    normal N(0.0, 1.0);

    double distance_to_default = (log(a0 / L) + (mu_a - 0.5 * pow(sigma_a, 2)) * t) / (sigma_a * sqrt(t));
    return cdf(N, -distance_to_default);
}

double get_min_capital_ROL(
    const double y,
    const double p,
    const double i
) {
    return ((1+y)*(1+p) - (1+i))/((1+i)*(1+y));
}

double get_returns_with_put(
    const double y,
    const double y_var,
    const double put,
    const double r
) {
    // Determine historic investment mean return
    const double sigma2_pt = log(1 + y_var/pow(y, 2));
    const double mu_pt = log(pow(y, 2) / sqrt(y_var + pow(y, 2)));

    // Determine the mean return with put protection
    normal norm(0.0, 1.0);
    const double zeta = (log(1+r) - mu_pt) / sqrt(sigma2_pt);
    double i_trunc = (1+r)*cdf(norm, zeta) + exp(mu_pt + sigma2_pt/2) - 1;

    return i_trunc;
}
