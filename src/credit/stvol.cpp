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
    auto realIntegrand = [this](double _phi) {
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
    auto realIntegrand = [this](double _phi) {
        auto numerator = charFn(_phi);
        auto denominator = 1i * _phi * pow(K, 1i * _phi);
        return real(numerator/denominator);
    };

    auto integrated = trapezoidal(realIntegrand, 1e-3, 1e3);
    return (0.5 + M_1_PI * integrated);
}

std::unique_ptr<StVol::Underlying> StVol::market_calibration(
    std::vector<HestonCallMdl>& hMdls, 
    const std::vector<double> market_prices,
    std::vector<double> initial_guess
)
{
    // Define a separate struct to simplify loading fixed data
    struct ExtData 
    {
        double S0;
        std::vector<double> rf;
        std::vector<double> strike;
        std::vector<double> t;
        std::vector<double> mdlPrices;
        std::vector<double> expPrices;
    };

    // NLopt requires objective functions to have following signature
    // Use the void* data part to pass in necessary data for expected and observed
    auto obj = [](const std::vector<double> &x, std::vector<double> &grad, void *data)
    {
        // Cast the void ptr to data
        auto* extData = static_cast<ExtData*>(data);
        const auto n_obs = extData->mdlPrices.size();
        auto err = 0.0;
        
        for (size_t j = 0; j < n_obs; ++j)
        {
            auto U = std::make_unique<StVol::Underlying>();

            // Fixed params
            U->S0 = extData->S0;

            // Optimization variables to pack into U-> members
            // Note: Structured bindings such as below are not permitted
            // [U->v0, U->alpha, U->vTheta, U->vSig, U->vLambda, U->rho] = x;
            // Instead use pointer-based loop to get around this
            std::vector<double*> U_elems = {&U->v0, &U->alpha, &U->vTheta, &U->vSig, &U->vLambda, &U->rho, &U->rf};

            for (size_t j = 0; j < U_elems.size(); ++j)
                *U_elems[j] = x[j];

            auto strike = extData->strike[j];
            auto t = extData->t[j];

            // Create Heston Model and evaluate the model price
            StVol::HestonCallMdl mdl(std::move(U), strike, t);
            mdl.calc_option_price();
            extData->mdlPrices[j] = mdl.get_option_price();

            err += pow(extData->mdlPrices[j] - extData->expPrices[j], 2) / n_obs;
        }

        return err;
    };

    // Prepare the data for optimization
    ExtData mdlData = {
        .S0 = hMdls[0].get_underlying().S0,
        .expPrices = market_prices
    };

    for (size_t j = 0; j < hMdls.size(); ++j)
    {
        mdlData.strike.push_back(hMdls[j].get_strike()) ;
        mdlData.t.push_back(hMdls[j].get_maturity());
        mdlData.rf.push_back(hMdls[j].get_underlying().rf);

        hMdls[j].calc_option_price();
        mdlData.mdlPrices.push_back(hMdls[j].get_option_price());
    }

    // Minimize the square error for calibration
    // Opt technique selected COBYLA (Constrained Opt. By Lin. Apprx)
    // Adv: no derivative required, good for non-convex smooth problems
    nlopt::opt optimizer(nlopt::LN_COBYLA, initial_guess.size()); 
    optimizer.set_min_objective(obj, &mdlData);
    optimizer.set_ftol_abs(1e-3); // Tolerance
    optimizer.set_maxeval(1000); // Maximum evaluation steps

    auto minf = 100.0;

    auto opt_code = optimizer.optimize(initial_guess, minf);
    std::cout << "Optimizer result code: " << opt_code << "\n\n";

    // Tidy the result by putting it into a new Underlying object
    auto result = std::make_unique<StVol::Underlying>();
    result->S0 = mdlData.S0;
    result->v0 = initial_guess[0];
    result->alpha = initial_guess[1];
    result->vTheta = initial_guess[2];
    result->vSig = initial_guess[3];
    result->vLambda = initial_guess[4];
    result->rho = initial_guess[5];
    result->rf = initial_guess[6];

    return result;
}
