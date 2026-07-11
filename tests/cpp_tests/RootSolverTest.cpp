// Testing file for the bracketing root solver (TOMS 748) - Boost Test Suite
#include "../../src/credit/utils.hpp"

#include <stdexcept>

// Includes Boost Test framework header. This macro handles dependencies and includes the runner
#define BOOST_TEST_MODULE RootSolverTestSuite
#include <boost/test/unit_test.hpp>

// Function definition remains the same
double f(double x) {
    return (3 * x * x + 2 * x - 8);
}

// Test case for root solver convergence verification
BOOST_AUTO_TEST_CASE(TestBracketRootConvergence) {
    // Detail original procedural steps for context/logging
    BOOST_TEST_MESSAGE("Solving 3x^2 +2x - 8 = 0 for x.");

    // Calculate required values
    const double f_exact = (-2 + std::sqrt(2 * 2 - 4 * 3 * (-8))) / (2 * 3);
    auto f_approx = bracket_root(f, 1.0);
    double error = (f_approx - f_exact) / f_exact;

    // Output results for comparison
    BOOST_TEST_MESSAGE("Exact solution     = " << f_exact);
    BOOST_TEST_MESSAGE("Solver solution    = " << f_approx);
    BOOST_TEST_MESSAGE("Percentage error   = " << error * 100 << "%");

    // Assertion based on the original successful pass condition
    BOOST_CHECK(std::abs(error) < 1e-3);
}

// A flat function has no sign change to bracket - the old secant iteration
// divided by zero here; the bracketing solver must instead fail loudly
BOOST_AUTO_TEST_CASE(TestFlatFunctionThrows) {
    BOOST_CHECK_THROW(bracket_root([](double) { return 1.0; }, 1.0), std::runtime_error);
}

// Root far outside the initial [x0/2, 2*x0] neighborhood - verifies the
// geometric bracket expansion
BOOST_AUTO_TEST_CASE(TestBracketExpansion) {
    auto root = bracket_root([](double x) { return x - 1000.0; }, 1.0);
    BOOST_CHECK(std::abs(root - 1000.0) < 1e-2);
}
