// Test implementation for function option_price

#ifndef TOL
#define TOL 1e-3
#endif

#include <iostream>
#include <math.h>
#include <vector>

#include "../source/risk_neutral.h"

double init = 100.0;
double strike = 120.0;
double interest = 0.05;
double volatility = 0.10;
double duration = 0.75;

double call_price, call_phi1, call_phi2;
double dividend = 0.001;  // For the put only
double put_price, put_phi1, put_phi2;

// Values calculated in Python - we expect these to be returned.
double expected_call_price = 0.18199750501184386;
double expected_call_Phi1 = 0.05166137589623915;
double expected_call_Phi2 = 0.04312161692483895;

double expected_put_price = 15.836460562018175;
double expected_put_Phi1 = 0.949248908843968;
double expected_put_Phi2 = 0.9576656456260665;

int main() {
    std::cout << "Testing function risk_neutral::vanilla_option_price." << std::endl;

    vanilla_option_price(init, strike, interest, volatility, duration, call_price, call_phi1, call_phi2);

    std::cout << "Call price    =    " << call_price << std::endl;
    std::cout << "Phi1          =    " << call_phi1 << std::endl;
    std::cout << "Phi2          =    " << call_phi2 << std::endl;

    if (abs(call_price - expected_call_price) <= TOL)
        std::cout << "Call price as expected" << std::endl;
    else
        std::cout << "Check call price!!\n";

    if (abs(call_phi1 - expected_call_Phi1) <= TOL)
        std::cout << "Phi1 as expected" << std::endl;
    else
        std::cout << "Check Phi1!!\n";

    if (abs(call_phi2 - expected_call_Phi2) <= TOL)
        std::cout << "Phi2 as expected" << std::endl;
    else
        std::cout << "Check Phi2!!\n";

    vanilla_option_price(init, strike, interest, volatility, duration, put_price, put_phi1, put_phi2, false, dividend);
    std::cout << "Put price     =   " << put_price << std::endl;
    std::cout << "Phi1          =    " << put_phi1 << std::endl;
    std::cout << "Phi2          =    " << put_phi2 << std::endl;

    if (abs(put_price - expected_put_price) <= TOL)
        std::cout << "Put price as expected" << std::endl;
    else
        std::cout << "Check put price!!\n";

    if (abs(put_phi1 - expected_put_Phi1) <= TOL)
        std::cout << "Phi1 as expected" << std::endl;
    else
        std::cout << "Check Phi1!!\n";

    if (abs(put_phi2 - expected_put_Phi2) <= TOL)
        std::cout << "Phi2 as expected" << std::endl;
    else
        std::cout << "Check Phi2!!\n";

    return 0;
}
