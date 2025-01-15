// Testing file for Heston model via theory
#include <cmath>
#include <iostream>
#include <memory>
#include <vector>
#include "../../src/credit/stvol.hpp"

using std::pow;
using std::sqrt;

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
    std::cout << "Volatility volatility: " << U->vSig << "\n";
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

bool test_Heston_calibration()
{
    std::vector<StVol::HestonCallMdl> mdls;
    std::vector<double> optnPrices;

    // Generate some artifical option prices
    std::cout << "Generating artificial option prices...\n";

    // Expected values
    auto S0 = 100.;
    auto v0 = 0.1;
    auto alpha = 1.5768;
    auto vTheta = 0.0398;
    auto vLambda = 0.575;
    auto vSig = 0.3;
    auto rho = -0.5711;

    // Set up initial values (expected values with some shifts)
    std::vector<double> x0 = {v0+0.05, alpha-0.01, vTheta+0.02, vSig+0.05, vLambda-0.01, rho+0.001, 0.03};

    for (int j = -10; j < 10; j++) 
    {
        auto dj = static_cast<double>(j);
        auto U = std::make_unique<StVol::Underlying>();
        U->S0 = S0;
        U->v0 = x0[0] * (1 + dj/100);
        U->alpha = x0[1];
        U->vTheta = x0[2] * (1 + 2*dj/100);
        U->vSig = x0[3] * (1 + dj/100);
        U->vLambda = x0[4] * (1 + 2*dj/100);
        U->rho = x0[5] * (1 + dj/100);
        U->rf = x0[6];
        auto strike = 100.;

        StVol::HestonCallMdl mdl(std::move(U), strike);
        mdl.calc_option_price();
        auto mdl_price = mdl.get_option_price();
        
        mdls.push_back(std::move(mdl));
        optnPrices.push_back(mdl_price);
    }
    std::cout << "Optimization step...\n";

    // Using the above, try to fit a Heston model to see if we recover params
    auto V = StVol::market_calibration(mdls, optnPrices, x0);
    std::cout << "Fitted result:\n";
    std::cout << "Spot price:            " << V->S0 << "\n";
    std::cout << "Spot volatility:       " << V->v0 << "\n";
    std::cout << "Risk-free rate:        " << V->rf << "\n";
    std::cout << "Mean-reversion rate:   " << V->alpha << "\n";
    std::cout << "Long-term mean var:    " << V->vTheta << "\n";
    std::cout << "Volatility volatility: " << V->vSig << "\n";
    std::cout << "Market price of vol:   " << V->vLambda << "\n";
    std::cout << "Correlation vol/price: " << V->rho << "\n";

    auto error2 = pow(V->v0 - v0, 2) / v0 + pow(V->alpha - alpha, 2) / alpha + pow(V->vTheta - vTheta, 2) / vTheta + pow(V->vLambda - vLambda, 2) / vLambda + pow(V->vSig - vSig, 2) / vSig + pow(V->rho - rho, 2) / rho;

    std::cout << "Avg percentage error:  " << sqrt(error2) * 100 / x0.size() << "%\n";

    StVol::HestonCallMdl mdl(std::move(V), 100.0);
    std::cout << "Fitted risk-neutral P2:" << mdl.get_rn_exercise_probability() << "\n";

    return (std::abs(sqrt(error2)) < 0.1 * x0.size());
}

int main() 
{
    auto testPass = true;

    std::cout << "Running test_Heston_pricing...\n\n";
    testPass &= test_Heston_pricing();
    std::cout << "\ntest_Heston_pricing result: " << (testPass?"Pass":"Fail") << "\n\n";

    // Final result
    std::cout << "Running test_Heston_calibration...\n\n";
    testPass &= test_Heston_calibration();
    std::cout << "\ntest_Heston_calibration result: " << (testPass?"Pass":"Fail") << "\n\n";
    std::cout << "All tests result: " << (testPass?"Pass":"Fail") << "\n";

    return (testPass?0:1);
}
