#include <cstdlib>
#ifndef TOL
#define TOL 1e-5
#endif

#include <iostream>
#include <cmath>
using std::cout;

#include "../source/template_utils.hpp"

auto f(double x) -> double {
    return (3 * x * x + 2 * x - 8);
}


int main() {
    cout << "Testing secant function by solving for one of the roots of equation:\n";
    cout << "3x^2 +2x - 8 = 0\n";
    auto f_exact = (-2 + std::sqrt(2*2 - 4 * 3 * (-8))) / (2 * 3);
    cout << "\n";
    cout << "Exact solution x1      = " << f_exact << "\n";

    auto f_approx = secant_root(f, 1.0);
    cout << "Secant solution xbar1  = " << f_approx << "\n";

    cout << "Result                 = " << ((std::abs(f_exact - f_approx) < TOL)?"Pass":"Fail") << "\n";

    return 0;
}
