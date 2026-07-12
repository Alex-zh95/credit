// Testing file for Heston model via theory - Adapted for Boost Test Suite
#include "../../src/credit/params.hpp"
#include "../../src/credit/stvol.hpp"
#include <cmath>
#include <iostream>
#include <memory>
#include <vector>

// Includes Boost Test framework header. This macro handles dependencies and includes the runner
// library.
#define BOOST_TEST_MODULE HestonTestSuite
#include <boost/test/unit_test.hpp>

// Test case 1: Heston pricing model validation
BOOST_AUTO_TEST_CASE(TestHestonPricing) {
    // Setup parameters based on original implementation
    StVol::HestonUnderlying U;
    U.S0 = 100.;
    U.v0 = 0.1;
    U.r = 0.03; // risk-free rate inherited from StandardUnderlying
    U.alpha = 1.5768;
    U.vTheta = 0.0398;
    U.vSig = 0.3;
    U.vLambda = 0.575;
    U.rho = -0.5711;

    // Expected outcome comparison - reference value computed independently
    // from the classic Heston (1993) two-probability formulation
    // C = S0*P1 - K*e^(-rt)*P2 with the same parameters
    const double expected_price = 10.820520;
    const double tolerance = 1e-3;

    // Running the pricing model
    auto strike = 100.;
    StVol::HestonCallMdl mdl(std::move(U), strike);
    mdl.calc_option_price();

    auto mdl_price = mdl.get_option_price();
    double error = (mdl_price - expected_price) / expected_price;

    // Use BOOST_CHECK for assertion
    BOOST_TEST_MESSAGE("Expected price: " << expected_price);
    BOOST_TEST_MESSAGE("Model price: " << mdl_price);
    BOOST_TEST_MESSAGE("Percentage error: " << error * 100 << "%");
    BOOST_TEST_MESSAGE("Risk-neutral P2: " << mdl.get_rn_exercise_probability());
    BOOST_TEST_MESSAGE("Delta (P1): " << mdl.get_delta());

    // Check if the error is within tolerance
    BOOST_CHECK(std::abs(error) < tolerance);

    // Delta equals P1 and the exercise probability equals P2 from the
    // two-probability formulation; reference values computed independently
    const double expected_delta = 0.638713;
    const double expected_p2 = 0.546664;
    BOOST_CHECK(std::abs(mdl.get_delta() - expected_delta) < tolerance);
    BOOST_CHECK(std::abs(mdl.get_rn_exercise_probability() - expected_p2) < tolerance);
}

// Test case 2: Heston asset volatility implied parameters calculation
BOOST_AUTO_TEST_CASE(TestHestonAssetVolatility) {
    // Alternative declaration: v0, alpha, vTheta, vSig, vLambda, rho
    std::vector<double> vols = {0.1, 1.5768, 0.0398, 0.3, 0.575, -0.5711};

    auto strike = 100.;
    auto rf = 0.03;

    StVol::HestonCallMdl mdl(StandardUnderlying(strike, strike, rf, 0.0), vols);

    auto V = StVol::HestonAssetVolatilityImplied(mdl, 150., 50., 1.5);

    // Outputting values for inspection (equivalent to original console output)
    BOOST_TEST_MESSAGE("Spot price: " << V.S0);
    BOOST_TEST_MESSAGE("Spot volatility: " << V.v0);
    // Add more checks/assertions here if specific values were expected in this scenario.
}

// Utility function for calibration demo (kept outside of main test flow)
void Heston_calibration_demo() {
    const auto S0 = 228.26;

    std::vector<double> Tau = {0.02, 0.02, 0.02, 0.04, 0.057};
    std::vector<double> Ks = {115, 120, 130, 115, 185};
    std::vector<double> Ps = {113.4, 108.4, 98.4, 114.6, 44.5};
    std::vector<double> Rfs = {0.0446, 0.0446, 0.0446, 0.0444, 0.0443};
    std::vector<double> trade_vols(Ps.size(), 1);

    auto U = StVol::fitHeston(S0, Ks, Rfs, Tau, Ps, trade_vols);

    // Outputting calibration results
    std::cout << "--- Heston Calibration Demo ---\n";
    std::cout << "Spot price:            " << U.S0 << "\n";
    std::cout << "Spot volatility:       " << U.v0 << "\n";
    std::cout << "Risk-free rate:        " << U.r << "\n";
    std::cout << "Mean-reversion rate:   " << U.alpha << "\n";
    std::cout << "Long-term mean var:    " << U.vTheta << "\n";
    std::cout << "Vol of vol:            " << U.vSig << "\n";
    std::cout << "Market price of vol:   " << U.vLambda << "\n";
    std::cout << "Correlation vol/price: " << U.rho << "\n";
}

// Note: The original int main() is replaced by the Boost automatic test runner setup.
// If this file needs to be runnable stand-alone without requiring CTest/Boost linkage,
// one might keep a minimal main(), but for proper CTest integration, the BOOST_TEST_MODULE is
// preferred. We rely on the Boost Test framework configuration for execution via CTest.
