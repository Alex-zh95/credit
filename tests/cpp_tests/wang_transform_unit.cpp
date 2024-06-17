#ifndef TOL
#define TOL 1e-3
#endif

#include <iostream>
#include <vector>

#include <boost/math/distributions.hpp>
using boost::math::normal;
#include "../../src/credit/risk_neutral.hpp"

template <class T>
int array_size(T& arr) {
    return (sizeof(arr) / sizeof(arr[0]));
}

int main() {
    // Define a distribution and discretize it
    normal N(5.1, 2.0);

    // Define some possible values
    double values[] = {2.3, 5.5, 9.1, -2.3, 2.0};
    const double sr = 1.2; // Sharpe ratio
    std::vector<double> P;

    double p;

    std::cout << "Testing function risk_neutral::wang_transform.\n";
    std::cout << "Taking certain quantiles of a normal distribution, we have:\n\t(";
    for (auto i = 0; i < array_size(values); ++i) {
        p = cdf(N, values[i]);
        P.insert(P.end(), p);

        std::cout << p << ((i == array_size(values) - 1)?")\n":", ");
    }

    std::cout << "Performing forward Wang transform...\n";
    std::vector<double> Q = wang_transform(P, sr);

    std::cout << "Output of transform:\t(";
    for (decltype(Q.size()) i = 0; i < Q.size(); ++i) {
        std::cout << Q[i] << ((i == Q.size() - 1)?")\n":", ");
    }

    std::cout << "Reverse the above to see if we recover the original...\n";
    std::vector<double> P_inv = wang_transform(Q, sr, true);

    std::cout << "Reverse:\t(";
    for (decltype(P_inv.size()) i = 0; i < P_inv.size(); ++i) {
        std::cout << P[i] << ((i == P_inv.size() - 1)?")\n":", ");
    }
}
