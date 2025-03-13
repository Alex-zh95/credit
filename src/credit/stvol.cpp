#include <cmath>
#include <complex>
#include <cstddef>
#include <memory>

using std::pow;
using std::exp;
using std::sqrt;
using std::real;
using namespace std::complex_literals;

#include <boost/math/quadrature/trapezoidal.hpp>
using boost::math::quadrature::trapezoidal;

#include <nlopt.hpp>

#include "stvol.hpp"
#include "template_utils.hpp"

std::complex<double> StVol::HestonCallMdl::charFn(std::complex<double> phi)
{
    auto a = underlying->alpha * underlying->vTheta;
    auto b = underlying->alpha + underlying->vLambda;

    // Parameter d given phi and b
    auto d = sqrt( pow(underlying->rho * underlying->vSig * phi * 1i - b, 2) + (phi * 1i + pow(phi, 2)) * pow(underlying->vSig, 2) );

    // Parameter g given g given phi, b and d
    auto g = (b - underlying->rho * underlying->vSig * phi * 1i + d) / (b - underlying->rho * underlying->vSig * phi * 1i - d);

    // Characteristic fn via the above components. This takes the form
    // exp(x1) * x2 * exp(x3)
    auto exp_x1 = exp(underlying->rf * phi * 1i * t);
    auto x2 = pow(underlying->S0, (phi * 1i)) * pow( (1. - g * exp(d * t)) / (1. - g), (-2. * a / pow(underlying->vSig, 2)) );
    auto exp_x2 = exp(
        a * t * (b - underlying->rho * underlying->vSig * phi * 1i + d) / pow(underlying->vSig, 2)
            + underlying->v0 * (b - underlying->rho * underlying->vSig * phi * 1i + d) * ( (1. - exp(d * t)) / (1. - g * exp(d * t)) )
            / (pow(underlying->vSig, 2))
    );

    return exp_x1 * x2 * exp_x2;
}

std::complex<double> StVol::HestonCallMdl::integrand(double phi)
{
    auto numerator = exp(underlying->rf * t) * charFn(phi - 1i) - K * charFn(phi);
    auto denominator = 1i * phi * pow(K, 1i * phi);

    return numerator/denominator;
}

void StVol::HestonCallMdl::calc_option_price()
{
    // Get the real part of the integrand
    auto realIntegrand = [this](double _phi) { return real(integrand(_phi)); };

    // Price is given by 0.5*(S0 - Ke^(-rf*t)) + 1/pi * integral(realIntegrand)
    auto integrated = trapezoidal(realIntegrand, 1e-3, 1e3);
    P = 0.5 * (underlying->S0 - K * exp(-underlying->rf * t)) + M_1_PI * integrated;
}

double StVol::HestonCallMdl::get_delta()
{
    auto realIntegrand = [this](double _phi) 
    {
        auto phiShift = _phi - 1i;
        auto numerator = charFn(phiShift);
        auto denominator = 1i * _phi * pow(K, 1i * _phi);
        return real(numerator/denominator);
    };

    auto integrated = trapezoidal(realIntegrand, 1e-3, 1e3);
    return (0.5 + M_1_PI * integrated);
}

double StVol::HestonCallMdl::get_rn_exercise_probability()
{
    auto realIntegrand = [this](double _phi) 
    {
        auto numerator = charFn(_phi);
        auto denominator = 1i * _phi * pow(K, 1i * _phi);
        return real(numerator/denominator);
    };

    auto integrated = trapezoidal(realIntegrand, 1e-3, 1e3);
    return (0.5 + M_1_PI * integrated);
}


std::unique_ptr<StVol::Underlying> StVol::fitHeston(double spot_price, std::vector<double> strikes, std::vector<double> r, std::vector<double> maturities, std::vector<double> market_prices, std::vector<double> trade_volumes)
{
    const auto nCalls = market_prices.size();

    // Extract the available volatility surface and set into parameters
    struct Params
    {
        std::vector<double> K;
        std::vector<double> t;
        std::vector<double> P;
        std::vector<double> rf;
        std::vector<double> volume;
        double S0;

    };

    Params params;
    
    params.K.reserve(nCalls);
    params.t.reserve(nCalls);
    params.P.reserve(nCalls);
    params.rf.reserve(nCalls);
    params.volume.reserve(nCalls);

    for (size_t i = 0; i < nCalls; ++i)
    {
        params.K.push_back(strikes.at(i));
        params.t.push_back(maturities.at(i));
        params.P.push_back(market_prices.at(i));
        params.rf.push_back(r.at(i));
        params.volume.push_back(trade_volumes.at(i));
    }

    params.S0 = spot_price;

    // Set up the variables to optimize in a vector, with the following order:
    // v0, alpha, vTheta, vSig, vLambda, rho
    std::vector<double> xVars = {0.1, 3.0, 0.05, 0.3, 0.03, -0.1};

    // Apply upper and lower bounds
    std::vector<double> xUb = {0.5, 5, 0.1, 1, 1, 1};
    std::vector<double> xLb = {1e-3, 1e-3, 1e-3, 1e-2, -1, -1};

    // NLopt requires objective functions to use following signature:
    // (const std::vector<double> &x, std::vector<double> &grad, void *data)
    // with x being the input vars to optimize, grad = gradient and data containing params
    auto square_err = [](const std::vector<double> &x, std::vector<double> &grad, void *data)
    {
        auto error = 0.0;

        // Cast void ptr to Params struct
        auto* parameters = static_cast<Params*>(data);

        // Iterate and count the errors
        const auto nOptions = parameters->P.size();
        auto total_volume = 0.0;
        for (size_t i = 0; i < nOptions; ++i)
        {
            auto curActualPrice = parameters->P.at(i);
            auto curStrike = parameters->K.at(i);
            auto curMaturity = parameters->t.at(i);
            auto curVolume = parameters->volume.at(i);

            auto U = std::make_unique<StVol::Underlying>();
            U->S0 = parameters->S0;
            U->v0 = x.at(0);
            U->alpha = x.at(1);
            U->vTheta = x.at(2);
            U->vSig = x.at(3);
            U->vLambda = x.at(4);
            U->rho = x.at(5);
            U->rf = parameters->rf.at(i);

            StVol::HestonCallMdl mdl(std::move(U), curStrike, curMaturity);
            mdl.calc_option_price();

            error += pow(curActualPrice - mdl.get_option_price(), 2) * curVolume;
            total_volume += curVolume;
        }

        return (error / total_volume);
    };

    // As we do not provide a gradient, require a deriv-free algo as no fd-approx implemented
    nlopt::opt optimizer(nlopt::LN_NELDERMEAD, xVars.size());
    optimizer.set_min_objective(square_err, &params);
    optimizer.set_xtol_abs(TOL);
    optimizer.set_maxeval(1e4);
    optimizer.set_upper_bounds(xUb);
    optimizer.set_lower_bounds(xLb);

    double minf;
    optimizer.optimize(xVars, minf);

    auto result = std::make_unique<StVol::Underlying>();
    result->S0 = spot_price;
    result->v0 = xVars.at(0);
    result->alpha = xVars.at(1);
    result->vTheta = xVars.at(2);
    result->vSig = xVars.at(3);
    result->vLambda = xVars.at(4);
    result->rho = xVars.at(5);
    result->rf = r.at(0);
    return result;
}


std::unique_ptr<StVol::Underlying> HestonAssetVolatilitySimulated(StVol::HestonCallMdl& mdl, double asset, double debt, double maturity)
{
    /* TODO:
     *
     * We want to implement a function that converts equity volatility into asset
     * volatility, in a similar vein to how the Merton model did this.
     *
     * We cannot use the direct relationship of sig_E * E = Delta * sig_A * A
     *
     * because volatility itself is stochastic here, as well as the E and A parts.
     *
     * We will make use of simulations here to do this, which will be more
     * computationally taxing but we use the power of C++ here. If possible we 
     * will also explore parallelization techniques but this might be difficult
     * if the process overall is iterative like it already is for Merton.
     *
     * Idea to implement here:
     *
     * - we use Heston to imply the equity value just as in Merton,
     * - from market data, we already have equity value, asset value and 
     *   the Heston model already calibrated to attain implied vol surface (equity),
     * 
     * Hence we will need as inputs:
     *
     * - a reference to a Heston Mdl object for equity
     * - information on asset (asset price (underlying), debt (strike))
     * 
     * The steps reqired therefore are:
     *
     * - forward solving of the Heston model to attain a price of equity, using asset price etc.
     * - compare this price of equity against the actual equity, i.e. mdl->get_underlying.S0
     *- minimize the square error here in a similar vein to HestonFit to parametrize the Heston parameters
     */

    // Variables to optimize in a vetor, with following order:
    // v0, alpha, vTheta, vSig, vLambda, rho
    std::vector<double> xVars = {
        mdl.get_underlying().v0,
        mdl.get_underlying().alpha,
        mdl.get_underlying().vTheta,
        mdl.get_underlying().vSig,
        mdl.get_underlying().vLambda,
        mdl.get_underlying().rho
    };

    // Apply upper and lower bounds
    std::vector<double> xUb = {0.5, 5, 0.1, 1, 1, 1};
    std::vector<double> xLb = {1e-3, 1e-3, 1e-3, 1e-2, -1, -1};

    // Parameters
    struct Params
    {
        double A0;
        double debt;
        double maturity;
        double actualEquity;
        double rf;
    };

    Params params;

    params.A0 = asset;
    params.debt = debt;
    params.maturity = maturity;
    params.actualEquity = mdl.get_underlying().S0;
    params.rf = mdl.get_underlying().rf;

    auto square_err = [](const std::vector<double> &x, std::vector<double> &grad, void* data)
    {
        auto* parameters = static_cast<Params*>(data);

        auto MertonUnderlying = std::make_unique<StVol::Underlying>();
        MertonUnderlying->S0 = parameters->A0;
        MertonUnderlying->v0 = x.at(0);
        MertonUnderlying->alpha = x.at(1);
        MertonUnderlying->vTheta = x.at(2);
        MertonUnderlying->vSig = x.at(3);
        MertonUnderlying->vLambda = x.at(4);
        MertonUnderlying->rho = x.at(5);
        MertonUnderlying->rf = parameters->rf;
        
        StVol::HestonCallMdl HestonMerton(std::move(MertonUnderlying), parameters->debt, parameters->maturity);

        HestonMerton.calc_option_price();
        auto calcEquity = HestonMerton.get_option_price();

        return pow(calcEquity - parameters->actualEquity, 2);
    };

    nlopt::opt optimizer(nlopt::LN_NELDERMEAD, xVars.size());
    optimizer.set_min_objective(square_err, &params);
    optimizer.set_xtol_abs(TOL);
    optimizer.set_maxeval(1e4);
    optimizer.set_upper_bounds(xUb);
    optimizer.set_lower_bounds(xLb);

    double minf;
    optimizer.optimize(xVars, minf);

    auto result = std::make_unique<StVol::Underlying>();
    result->S0 = asset;             // Asset price
    result->v0 = xVars.at(0);       // Spot asset volatility
    result->alpha = xVars.at(1);    // Mean reversion rate
    result->vTheta = xVars.at(2);   // Long-term asset volatility
    result->vSig = xVars.at(3);     // Volatility of asset volatility
    result->vLambda = xVars.at(4);  // Market price of asset volatility
    result->rho = xVars.at(5);      // Asset price to volatility correlation
    result->rf = params.rf;         // Risk-free rate
    return result;
}
