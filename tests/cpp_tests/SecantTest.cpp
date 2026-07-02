// Testing file for Secant method implementation - Adapted for Boost Test Suite
#include "../../src/credit/utils.hpp"

// Includes Boost Test framework header. This macro handles dependencies and includes the runner
#define BOOST_TEST_MODULE SecantTestSuite
#include <boost/test/unit_test.hpp>

// Function definition remains the same
double f(double x) {
    return (3 * x * x + 2 * x - 8);
}

// Test case for Secant method convergence verification
BOOST_AUTO_TEST_CASE(TestSecantConvergence) {
    // Detail original procedural steps for context/logging
    BOOST_TEST_MESSAGE("Solving 3x^2 +2x - 8 = 0 for x.");

    // Calculate required values
    const double f_exact = (-2 + std::sqrt(2 * 2 - 4 * 3 * (-8))) / (2 * 3);
    auto f_approx = secant_root(f, 1.0);
    double error = (f_approx - f_exact) / f_exact;

    // Output results for comparison
    BOOST_TEST_MESSAGE("Exact solution     = " << f_exact);
    BOOST_TEST_MESSAGE("Secant solution    = " << f_approx);
    BOOST_TEST_MESSAGE("Percentage error   = " << error * 100 << "%");

    // Assertion based on the original successful pass condition
    BOOST_CHECK(std::abs(error) < 1e-3);
}