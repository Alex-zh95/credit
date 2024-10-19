#include <iostream>
#include <cmath>

#include "../../src/credit/template_utils.hpp"

double f(double x) {
    return (3 * x * x + 2 * x - 8);
}

bool test_secant()
{
    std::cout << "Solving 3x^2 +2x - 8 = 0 for x.\n";

    auto f_exact = (-2 + std::sqrt(2*2 - 4 * 3 * (-8))) / (2 * 3);
    auto f_approx = secant_root(f, 1.0);
    auto error = (f_approx - f_exact) / f_exact;
    std::cout << "Exact solution     = " << f_exact << "\n";
    std::cout << "Secant solution    = " << f_approx << "\n";
    std::cout << "Percentage error   = " << error * 100 << "%\n";
    
    return (std::abs(error) < 1e-3);
}


int main() {
    auto testPass = true;

    std::cout << " Running test_secant...\n\n";
    testPass &= test_secant();
    std::cout << "\ntest_secant result: " << (testPass?"Pass":"Fail") << "\n\n";

    // Final result
    std::cout << "All tests result: " << (testPass?"Pass":"Fail") << "\n";

    return (testPass);
}
