#include <complex>

using std::log;
using std::pow;
using std::exp;
using std::sqrt;
using std::real;
using namespace std::complex_literals;

#include <boost/math/quadrature/trapezoidal.hpp>
using boost::math::quadrature::trapezoidal;

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
