// Implementation of prototypes in risk_neutral.h
#include <math.h>
#include <boost/math/distributions/normal.hpp>
using boost::math::normal;

#include "risk_neutral.h"

void option_price(
        double S0,
        double K,
        double r,
        double sigma,
        double t,
        double& price,
        double& Phi1,
        double& Phi2,
        bool call = true,
        double q = 0,
        ) {
    double d1, d2;
    double Phi1, Phi2;

    // Define N(0,1) object
    normal N(0.0, 1.0);

    d1 = (log(S0/K) + (r-q+0.5*pow(sigma,2.0) * t) / (sigma * sqrt(t)));
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
        double sharpe_ratio,
        bool inverse = false
        ) {
    // Create a copy of the input and manipulate this instead
    std::vector<double> Q = P;

    // Define N(0,1) object
    normal N(0.0, 1.0);

    if (!inverse) {
        for(auto& p : Q) {
            p = cdf(N, quantile(N, p) + sharpe_ratio)
        }
    } else {
        for(auto& p : Q) {
            p = cdf(N, quantile(N, p) - sharpe_ratio)
        }
    }

    return Q;
}
