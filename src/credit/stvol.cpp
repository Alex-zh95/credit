#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <future>
#include <numbers>
#include <numeric>
#include <thread>
#include <vector>

using std::exp;
using std::pow;
using std::real;
using std::sqrt;
using namespace std::complex_literals;

#include <boost/math/quadrature/trapezoidal.hpp>
using boost::math::quadrature::trapezoidal;

#include <nlopt.hpp>

#include "risk_neutral.hpp"
#include "stvol.hpp"
#include "thread_pool.hpp"
#include "utils.hpp"

std::complex<double> StVol::HestonCallMdl::charFn(std::complex<double> phi) const {
    const auto t = underlying.t;

    auto a = underlying.alpha * underlying.vTheta;
    auto b = underlying.alpha + underlying.vLambda;

    // Parameter d given phi and b
    auto d = sqrt(pow(underlying.rho * underlying.vSig * phi * 1i - b, 2) +
                  (phi * 1i + pow(phi, 2)) * pow(underlying.vSig, 2));

    // Parameter g given g given phi, b and d
    auto g = (b - underlying.rho * underlying.vSig * phi * 1i + d) /
             (b - underlying.rho * underlying.vSig * phi * 1i - d);

    // Characteristic fn via the above components. This takes the form
    // exp(x1) * x2 * exp(x3)
    auto exp_x1 = exp(underlying.r * phi * 1i * t);
    auto x2 = pow(underlying.S0, (phi * 1i)) *
              pow((1. - g * exp(d * t)) / (1. - g), (-2. * a / pow(underlying.vSig, 2)));
    auto exp_x2 = exp(a * t * (b - underlying.rho * underlying.vSig * phi * 1i + d) /
                          pow(underlying.vSig, 2) +
                      underlying.v0 * (b - underlying.rho * underlying.vSig * phi * 1i + d) *
                          ((1. - exp(d * t)) / (1. - g * exp(d * t))) / (pow(underlying.vSig, 2)));

    return exp_x1 * x2 * exp_x2;
}

std::complex<double> StVol::HestonCallMdl::integrand(double phi) const {
    const auto K = underlying.K;

    auto numerator = exp(underlying.r * underlying.t) * charFn(phi - 1i) - K * charFn(phi);
    auto denominator = 1i * phi * pow(K, 1i * phi);

    return numerator / denominator;
}

void StVol::HestonCallMdl::calc_option_price() {
    // Get the real part of the integrand
    auto realIntegrand = [this](double _phi) { return real(integrand(_phi)); };

    // Price is given by 0.5*(S0 - Ke^(-rf*t)) + 1/pi * integral(realIntegrand)
    auto integrated = trapezoidal(realIntegrand, 1e-3, 1e3);
    P = 0.5 * (underlying.S0 - underlying.K * exp(-underlying.r * underlying.t)) +
        std::numbers::inv_pi * integrated;
}

double StVol::HestonCallMdl::get_delta() const {
    auto realIntegrand = [this](double _phi) {
        auto phiShift = _phi - 1i;
        auto numerator = charFn(phiShift);
        auto denominator = 1i * _phi * pow(underlying.K, 1i * _phi);
        return real(numerator / denominator);
    };

    auto integrated = trapezoidal(realIntegrand, 1e-3, 1e3);
    return (0.5 + std::numbers::inv_pi * integrated);
}

double StVol::HestonCallMdl::get_rn_exercise_probability() const {
    auto realIntegrand = [this](double _phi) {
        auto numerator = charFn(_phi);
        auto denominator = 1i * _phi * pow(underlying.K, 1i * _phi);
        return real(numerator / denominator);
    };

    auto integrated = trapezoidal(realIntegrand, 1e-3, 1e3);
    return (0.5 + std::numbers::inv_pi * integrated);
}

StVol::HestonUnderlying StVol::fitHeston(double spot_price, std::vector<double> strikes,
                                         std::vector<double> r, std::vector<double> maturities,
                                         std::vector<double> market_prices,
                                         std::vector<double> trade_volumes) {
    const auto nCalls = market_prices.size();

    // Extract the available volatility surface and set into parameters
    struct Params {
        std::vector<double> K;
        std::vector<double> t;
        std::vector<double> P;
        std::vector<double> rf;
        std::vector<double> volume;
        double S0;
        ThreadPool* pool;
        std::size_t nWorkers;
    };

    Params params;

    params.K.reserve(nCalls);
    params.t.reserve(nCalls);
    params.P.reserve(nCalls);
    params.rf.reserve(nCalls);
    params.volume.reserve(nCalls);

    for (size_t i = 0; i < nCalls; ++i) {
        params.K.push_back(strikes.at(i));
        params.t.push_back(maturities.at(i));
        params.P.push_back(market_prices.at(i));
        params.rf.push_back(r.at(i));
        params.volume.push_back(trade_volumes.at(i));
    }

    params.S0 = spot_price;

    // Workers are created once here and reused across every objective
    // evaluation - NLopt calls the objective thousands of times, so spawning
    // threads inside it (the previous design) paid creation/teardown each call.
    // No point running more workers than there are options to price.
    std::size_t nWorkers = std::thread::hardware_concurrency();
    if (nWorkers == 0) nWorkers = 1;
    nWorkers = std::min<std::size_t>(nWorkers, nCalls > 0 ? nCalls : 1);

    ThreadPool pool(nWorkers);
    params.pool = &pool;
    params.nWorkers = nWorkers;

    // Set up the variables to optimize in a vector, with the following order:
    // v0, alpha, vTheta, vSig, vLambda, rho
    std::vector<double> xVars = {0.1, 3.0, 0.05, 0.3, 0.03, -0.1};

    // Apply upper and lower bounds
    std::vector<double> xUb = {0.5, 5, 0.1, 1, 1, 1};
    std::vector<double> xLb = {1e-3, 1e-3, 1e-3, 1e-2, -1, -1};

    // NLopt requires objective functions to use following signature:
    // (const std::vector<double> &x, std::vector<double> &grad, void *data)
    // with x being the input vars to optimize, grad = gradient and data containing params
    auto square_err = [](const std::vector<double>& x, std::vector<double>& grad, void* data) {
        // Cast void ptr to Params struct
        auto* parameters = static_cast<Params*>(data);

        // Iterate and count the errors
        const auto nOptions = parameters->P.size();
        const auto nWorkers = parameters->nWorkers;
        const std::size_t chunkSize = nOptions / nWorkers;
        const std::size_t chunkRemainder = nOptions % nWorkers;

        std::vector<double> partialErr(nWorkers, 0.0);
        std::vector<double> partialVol(nWorkers, 0.0);

        // Dispatch one chunk per worker; the futures both signal completion and
        // surface any exception thrown while pricing (raw threads would terminate)
        std::vector<std::future<void>> pending;
        pending.reserve(nWorkers);

        for (std::size_t j = 0; j < nWorkers; ++j) {
            std::size_t start = j * chunkSize;
            std::size_t end = start + chunkSize + ((j == nWorkers - 1) ? (chunkRemainder) : 0);

            pending.push_back(parameters->pool->submit(
                [parameters, &x, &partialErr, &partialVol, j, start, end] {
                    for (std::size_t i = start; i < end; ++i) {
                        auto curActualPrice = parameters->P.at(i);
                        auto curStrike = parameters->K.at(i);
                        auto curMaturity = parameters->t.at(i);
                        auto curVolume = parameters->volume.at(i);
                        auto curRiskFree = parameters->rf.at(i);

                        auto mdl =
                            StVol::HestonCallMdl(StandardUnderlying(parameters->S0, curStrike,
                                                                    curRiskFree, 0.0, curMaturity),
                                                 x);
                        mdl.calc_option_price();

                        partialErr[j] +=
                            pow(curActualPrice - mdl.get_option_price(), 2) * curVolume;
                        partialVol[j] += curVolume;
                    }
                }));
        }

        for (auto& fut : pending) fut.get();

        auto error = std::reduce(partialErr.begin(), partialErr.end());
        auto total_volume = std::reduce(partialVol.begin(), partialVol.end());

        return (error / total_volume);
    };

    // As we do not provide a gradient, require a deriv-free algo as no fd-approx implemented
    nlopt::opt optimizer(nlopt::LN_NELDERMEAD, xVars.size());
    optimizer.set_min_objective(square_err, &params);
    optimizer.set_xtol_abs(TOL);
    optimizer.set_maxeval(10000);
    optimizer.set_upper_bounds(xUb);
    optimizer.set_lower_bounds(xLb);

    double minf;
    optimizer.optimize(xVars, minf);

    // xVars order (v0, alpha, vTheta, vSig, vLambda, rho) matches the
    // HestonUnderlying volatility-parameter constructor
    return StVol::HestonUnderlying(StandardUnderlying(spot_price, 0.0, r.at(0), 0.0), xVars);
}

StVol::HestonUnderlying StVol::HestonAssetVolatilityImplied(const StVol::HestonCallMdl& mdl,
                                                            double asset, double debt,
                                                            double maturity) {
    StVol::HestonUnderlying U;

    // Use the Ito-derived relationship sig_E * E = Delta * sig_A * A,
    // to derive the spot and long-term volatility
    U.v0 = get_asset_volatility(asset, mdl.get_underlying().v0, debt, mdl.get_underlying().r,
                                maturity);

    U.vTheta = get_asset_volatility(asset, mdl.get_underlying().vTheta, debt,
                                    mdl.get_underlying().r, maturity);

    // Keeping it simple, the relationship above states sig_E and sig_A
    // are linear proportional at time t. So for volatility of volatility,
    // use a simple scaling
    U.vSig = mdl.get_underlying().vSig * U.v0 / mdl.get_underlying().v0;

    // For correlation, we need the implied equity to asset delta
    // Simplify by using standard Black-Scholes
    auto rf = mdl.get_underlying().r;
    auto [implied_E, implied_asset_equity_delta, Phi2] =
        vanilla_option_price(StandardUnderlying(asset, debt, rf, rf, maturity));

    // Asset volatility to value correlation motivation:
    // If there is no leverage, then this equates that of equity volatility
    // Otherwise, the implied asset equity delta is the measure of change
    // in asset price against that of equity, so this can be used to dampen rho
    U.rho = mdl.get_underlying().rho / implied_asset_equity_delta;

    // Mean reversion assumed linear as for volatility of volatility
    U.alpha = mdl.get_underlying().alpha * U.v0 / mdl.get_underlying().v0;

    // For market price of volatility, as the asset volatility is unobservable,
    // and that the market price for the equity volatility can be distilled
    // via Sharpe ratio, assume that the market price for asset volatility
    // follows the same exchange rate as that for vSig or mean reversion
    U.vLambda = mdl.get_underlying().vLambda * U.v0 / mdl.get_underlying().v0;

    // Remaining params (inherited from StandardUnderlying)
    U.S0 = asset;
    U.r = mdl.get_underlying().r;

    return U;
}
