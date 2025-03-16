// Implementation of prototypes in risk_neutral.h
#include <cmath>
using std::log;
using std::pow;
using std::exp;
using std::sqrt;

#include <boost/math/distributions/normal.hpp>
using boost::math::normal;

#include "utils.hpp"
#include "risk_neutral.hpp"

std::tuple<double, double, double> vanilla_option_price(
    const double S0,
    const double K,
    const double r,
    const double sigma,
    const double t,
    const bool call,
    const double q)
{
    auto Phi1 = 0.0;
    auto Phi2 = 0.0;
    auto price = 0.0;

    // Define N(0,1) object
    normal N(0.0, 1.0);

    auto d1 = (log(S0 / K) + (r - q + 0.5 * pow(sigma, 2.0)) * t) / (sigma * sqrt(t));
    auto d2 = d1 - sigma * sqrt(t);

    if (call)
    {
        Phi1 = cdf(N, d1);
        Phi2 = cdf(N, d2);

        price = S0 * exp(-q * t) * Phi1 - K * exp(-r * t) * Phi2;
    }
    else
    {
        Phi1 = cdf(N, -d1);
        Phi2 = cdf(N, -d2);

        price = K * exp(-r * t) * Phi2 - S0 * exp(-q * t) * Phi1;
    }

    return {price, Phi1, Phi2};
}

std::tuple<double, double> fpt_call_price(
    const double S0,
    const double K,
    const double r,
    const double sigma,
    const double t,
    const double gamma,
    const double q)
{
    normal N(0.0, 1.0);

    if (S0 < K)
        return {0.0, 0.0};

    auto Kt = K * exp(gamma * t);
    auto lambda = (r - q + 0.5 * pow(sigma, 2)) / pow(sigma, 2);
    auto y = log(Kt / S0) / (sigma * sqrt(t)) + lambda * sigma * sqrt(t);
    auto x = log(S0 / Kt) / (sigma * sqrt(t)) + lambda * sigma * sqrt(t);

    auto c_do = S0 * cdf(N, x) * exp(-q * t) - Kt * cdf(N, x - sigma * sqrt(t)) - S0 * exp(-q * t) * pow(Kt / S0, 2 * lambda) * cdf(N, y) + Kt * pow(Kt / S0, 2 * lambda - 2) * cdf(N, y - sigma * sqrt(t));

    // Exploiting the fact that a standard call can be decomposed
    // into sum of down-and-out and down-and-in calls
    auto [c, Phi1, Phi2] = vanilla_option_price(S0, K, r, sigma, t, true, q);
    auto c_di = c - c_do;

    return {c_do, c_di};
}

std::vector<double> wang_transform(
    std::vector<double> P,
    const double sharpe_ratio,
    const bool inverse)
{
    // Define N(0,1) object
    normal N(0.0, 1.0);

    for (auto &p : P)
        p = cdf(N, quantile(N, p) + (inverse ? -1 : 1) * sharpe_ratio);

    return P;
}

double get_asset_volatility(
    const double E,
    const double sigma_e,
    const double L,
    const double r,
    const double t,
    const int n_iter)
{
    // Check the trivial case: if no leverage, asset and equity volatility equate
    if (L < TOL)
        return sigma_e;

    // Initial guesses
    auto cur_asset = E;
    auto sigma_a = 0.5;

    // Iterative steps
    for (auto run_iter = 0; run_iter < n_iter; ++run_iter)
    {
        // Save down a previous result for sigma_a
        auto prev_sigma_a = sigma_a;

        // Update sigma_a with current asset value vector
        // From Itô Lemma, we have the following relationship: sigma_e * Equity / Asset = Phi1 * sigma_a
        auto [eq, Phi1, Phi2] = vanilla_option_price(cur_asset, L, r, sigma_a, t);
        sigma_a = sigma_e * E / cur_asset / Phi1;

        // Early stop if difference between previous and current sigma_iterations is below tolerance
        if (abs(sigma_a - prev_sigma_a) < TOL)
            break;

        // Imply current asset value using current guess of sigma_a
        auto fn = [&L, &r, &sigma_a, &t, &E](double _a)
        {
            auto [impl_equity, Phi1, Phi2] = vanilla_option_price(_a, L, r, sigma_a, t);
            return (E - impl_equity);
        };

        cur_asset = secant_root(fn, E);
    }

    return sigma_a;
}

double get_vanilla_default_probability(
    const double a0,
    const double rf,
    const double sigma_a,
    const double L,
    const double t)
{
    normal N(0.0, 1.0);

    auto distance_to_default = (log(a0 / L) + (rf - 0.5 * pow(sigma_a, 2)) * t) / (sigma_a * sqrt(t));
    return cdf(N, -distance_to_default);
}

double get_fpt_default_probability(
    const double a0,
    const double rf,
    const double sigma_a,
    const double L,
    const double q,
    const double gamma,
    const double t)
{
    normal N(0.0, 1.0);

    auto Lt = L * exp(gamma * t); // Default boundary via gamma growth rate at expiry

    auto d1 = (log(a0) - log(Lt) + (rf - q - 0.5 * pow(sigma_a, 2)) * t) / (sigma_a * sqrt(t));
    auto d2 = (-log(a0) - log(Lt) + (rf - q - 0.5 * pow(sigma_a, 2)) * t) / (sigma_a * sqrt(t));

    return 1 - (cdf(N, d1) - pow(a0 / Lt, 1 - 2 * (rf - q - gamma) / pow(sigma_a, 2)) * cdf(N, d2));
}

double get_min_capital_ROL(
    const double y,
    const double p,
    const double i)
{
    return ((1 + y) * (1 + p) - (1 + i)) / ((1 + i) * (1 + y));
}

double get_returns_with_put(
    const double y,
    const double y_var,
    const double put,
    const double r)
{
    // Determine historic investment mean return
    const auto sigma2_pt = log(1 + y_var / pow(y, 2));
    const auto mu_pt = log(pow(y, 2) / sqrt(y_var + pow(y, 2)));

    // Determine the mean return with put protection
    normal norm(0.0, 1.0);
    const auto zeta = (log(1 + r) - mu_pt) / sqrt(sigma2_pt);
    auto i_trunc = (1 + r) * cdf(norm, zeta) + exp(mu_pt + sigma2_pt / 2) - 1;

    return i_trunc;
}

double put_call_parity(    
    const double in_price,
    const double rf,
    const double S0,
    const double K,
    const double t,
    bool in_call)
{
    auto result = 0.0;

    if (in_call) 
        // Incoming price is call so return put
        result = in_price - S0 + K * exp(-rf * t);
    else 
        // Incoming price is put so return call
        result = S0 - K * exp(-rf * t) + in_price;

    return result;
}
