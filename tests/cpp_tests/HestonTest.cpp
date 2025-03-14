// Testing file for Heston model via theory
#include <cmath>
#include <iostream>
#include <memory>
#include <vector>
#include "../../src/credit/stvol.hpp"

bool test_Heston_pricing()
{
    auto U = std::make_unique<StVol::Underlying>();
    U->S0 = 100.;
    U->v0 = 0.1;
    U->rf = 0.03;
    U->alpha = 1.5768;
    U->vTheta = 0.0398;
    U->vSig = 0.3;
    U->vLambda = 0.575;
    U->rho = -0.5711;

    const auto expected_price = 11.540361819355377;
    
    std::cout << "Spot price:            " << U->S0 << "\n";
    std::cout << "Spot volatility:       " << U->v0 << "\n";
    std::cout << "Risk-free rate:        " << U->rf << "\n";
    std::cout << "Mean-reversion rate:   " << U->alpha << "\n";
    std::cout << "Long-term mean var:    " << U->vTheta << "\n";
    std::cout << "Vol of vol:            " << U->vSig << "\n";
    std::cout << "Market price of vol:   " << U->vLambda << "\n";
    std::cout << "Correlation vol/price: " << U->rho << "\n";

    auto strike = 100.;
    std::cout << "Strike price for call: " << strike << "\n\n";

    StVol::HestonCallMdl mdl(std::move(U), strike);
    mdl.calc_option_price();

    auto mdl_price = mdl.get_option_price();
    auto error = (mdl_price - expected_price) / expected_price;
    std::cout << "Expected price:        " << expected_price << "\n";
    std::cout << "Model price:           " << mdl_price << "\n";
    std::cout << "Percentage error:      " << error * 100 << "%\n";
    std::cout << "Risk-neutral P2:       " << mdl.get_rn_exercise_probability() << "\n";

    return (std::abs(error) < 1e-3);
}

void test_Heston_asset_vol()
{
    auto U = std::make_unique<StVol::Underlying>();
    U->S0 = 100.;
    U->v0 = 0.1;
    U->rf = 0.03;
    U->alpha = 1.5768;
    U->vTheta = 0.0398;
    U->vSig = 0.3;
    U->vLambda = 0.575;
    U->rho = -0.5711;

    auto strike = 100.;
    
    StVol::HestonCallMdl mdl(std::move(U), strike);

    auto V = StVol::HestonAssetVolatilityImplied(mdl, 150., 50., 1.5);

    std::cout << "Spot price:            " << V->S0 << "\n";
    std::cout << "Spot volatility:       " << V->v0 << "\n";
    std::cout << "Risk-free rate:        " << V->rf << "\n";
    std::cout << "Mean-reversion rate:   " << V->alpha << "\n";
    std::cout << "Long-term mean var:    " << V->vTheta << "\n";
    std::cout << "Volatility volatility: " << V->vSig << "\n";
    std::cout << "Market price of vol:   " << V->vLambda << "\n";
    std::cout << "Correlation vol/price: " << V->rho << "\n";
}

void basic_Heston_calibration_test()
{
    const auto S0 = 228.26;

    std::vector<double> Tau = {0.02, 0.02, 0.02, 0.04, 0.057};
    std::vector<double> Ks = {115, 120, 130, 115, 185};
    std::vector<double> Ps = {113.4, 108.4, 98.4, 114.6, 44.5};
    std::vector<double> Rfs = {0.0446, 0.0446, 0.0446, 0.0444, 0.0443};
    std::vector<double> trade_vols(Ps.size(), 1);       // Take the number of elements for Ps and write 1 in each

    auto U = StVol::fitHeston(S0, Ks, Rfs, Tau, Ps, trade_vols);
    std::cout << "Spot price:            " << U->S0 << "\n";
    std::cout << "Spot volatility:       " << U->v0 << "\n";
    std::cout << "Risk-free rate:        " << U->rf << "\n";
    std::cout << "Mean-reversion rate:   " << U->alpha << "\n";
    std::cout << "Long-term mean var:    " << U->vTheta << "\n";
    std::cout << "Volatility volatility: " << U->vSig << "\n";
    std::cout << "Market price of vol:   " << U->vLambda << "\n";
    std::cout << "Correlation vol/price: " << U->rho << "\n";
}

int main() 
{
    auto testPass = true;

    std::cout << "Running test_Heston_pricing...\n\n";
    testPass &= test_Heston_pricing();
    std::cout << "\ntest_Heston_pricing result: " << (testPass?"Pass":"Fail") << "\n\n";

    // std::cout << "All tests result: " << (testPass?"Pass":"Fail") << "\n";
    std::cout << "Testing calibration...\n";
    basic_Heston_calibration_test();

    // Testing Heston asset volatility
    std::cout << "\nTesting implied asset volatility...\n";
    test_Heston_asset_vol();
    return (testPass?0:1);
}
