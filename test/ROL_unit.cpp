// Test implementation of minimum ROL
#include <iostream>
using std::cout;

#include "../source/risk_neutral.h"

int main() {
    // *** Part 1 ***: Test with put option from market
    // Yahoo Finance (MSFT at 2024-02-16)
    const double target_growth = 0.12;  // Growth rate from Yahoo Finance
    const double put_strike = 0.05; // Strike at chosen risk-free rate

    // Corresponding put option rate from market corresponding to above
    const double put_price = 29.90 / 406.56;
    const double growth_variance = 0.2030; // Implied volatility from chosen option 

    cout << "Testing ROL calculations...with put option.\n";
    cout << "Put strike:            " << put_strike*100 << "%\n";
    cout << "Put price:             " << put_price*100 << "%\n";
    cout << "Target growth:         " << target_growth * 100 << "%\n";

    double i_trunc = get_returns_with_put(target_growth, growth_variance, put_price, put_strike);
    cout << "Growth with put:       " << i_trunc*100 << "%\n";
    double rol = get_min_ROL(target_growth, put_price, i_trunc);
    cout << "ROL = " << rol*100 << "%\n\n";

    // *** Part 2 ***: Test with no put option
    cout << "Testing ROL calculations...without pt option.\n";
    cout << "Using same target growth and risk-free rate.\n";
    double rol2 = get_min_ROL(target_growth, 0.0, put_strike);
    cout << "ROL = " << rol2*100 << "%\n";

    cout << "Option ROL < Swap ROL? " << ((rol <= rol)?"Yes":"No") << "\n\n";
    return 0;
}
