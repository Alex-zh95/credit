#include <cmath>
#include <complex>
#include <cstddef>
#include <memory>
#include <iostream>

using std::pow;
using std::exp;
using std::sqrt;
using std::real;
using namespace std::complex_literals;

#include <boost/math/quadrature/trapezoidal.hpp>
using boost::math::quadrature::trapezoidal;

#include <nlopt.hpp>

#include "stvol.hpp"

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


std::unique_ptr<StVol::Underlying> StVol::fitHeston(double spot_price, std::vector<double> strikes, std::vector<double> r, std::vector<double> maturities, std::vector<double> market_prices)
{
    const auto nCalls = market_prices.size();

    // Extract the available volatility surface and set into parameters
    struct Params
    {
        std::vector<double> K;
        std::vector<double> t;
        std::vector<double> P;
        std::vector<double> rf;
        double S0;

    };

    Params params;
    
    params.K.reserve(nCalls);
    params.t.reserve(nCalls);
    params.P.reserve(nCalls);
    params.rf.reserve(nCalls);

    for (size_t i = 0; i < nCalls; ++i)
    {
        params.K.push_back(strikes.at(i));
        params.t.push_back(maturities.at(i));
        params.P.push_back(market_prices.at(i));
        params.rf.push_back(r.at(i));
    }

    params.S0 = spot_price;

    // Set up the variables to optimize in a vector, with the following order:
    // v0, alpha, vTheta, vSig, vLambda, rho
    std::vector<double> xVars = {0.1, 3.0, 0.05, 0.3, 0.03, -0.1};

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
        for (size_t i = 0; i < nOptions; ++i)
        {
            auto curActualPrice = parameters->P.at(i);
            auto curStrike = parameters->K.at(i);
            auto curMaturity = parameters->t.at(i);

            auto U = std::make_unique<StVol::Underlying>();
            U->S0 = parameters->S0;
            U->v0 = x[0];
            U->alpha = x[1];
            U->vTheta = x[2];
            U->vSig = x[3];
            U->vLambda = x[4];
            U->rho = x[5];
            U->rf = parameters->rf.at(i);

            StVol::HestonCallMdl mdl(std::move(U), curStrike, curMaturity);
            mdl.calc_option_price();

            error += pow(curActualPrice - mdl.get_option_price(), 2);
        }

        return (error / nOptions);
    };

    nlopt::opt optimizer(nlopt::LD_SLSQP, xVars.size());
    optimizer.set_min_objective(square_err, &params);
    optimizer.set_xtol_abs(1e-3);
    optimizer.set_maxeval(1e4);

    optimizer.optimize(xVars);

    auto result = std::make_unique<StVol::Underlying>();
    result->S0 = spot_price;
    result->v0 = xVars[0];
    result->alpha = xVars[1];
    result->vTheta = xVars[2];
    result->vSig = xVars[3];
    result->vLambda = xVars[4];
    result->rho = xVars[5];
    result->rf = r.at(0);
    return result;
}
