// Testing file for Heston model via theory
#include <cmath>
#include <iostream>
#include <memory>
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

    return (std::abs(error) < 1e-3);
}

int main() 
{
    auto testPass = true;

    std::cout << "Running test_Heston_pricing...\n\n";
    testPass &= test_Heston_pricing();
    std::cout << "\ntest_Heston_pricing result: " << (testPass?"Pass":"Fail") << "\n\n";

    // Final result
    std::cout << "All tests result: " << (testPass?"Pass":"Fail") << "\n";

    return (testPass);
}
