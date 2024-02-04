// Test implementation for function option_price

#include <iostream>
#include "../source/risk_neutral.h"

double init = 100.0;
double strike = 120.0;
double interest = 0.05;
double volatility = 0.10;
double duration = 0.75;

double call_price, call_phi1, call_phi2;

int main() {
    std::cout << "Testing function risk_neutral::option_price." << std::endl;

    option_price(init, strike, interest, volatility, duration, call_price, call_phi1, call_phi2);

    std::cout << "Option price = " << call_price << std::endl;

    return 0;
}
