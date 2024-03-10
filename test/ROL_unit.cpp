// Test implementation of minimum ROL
#include <cstdlib>
#include <iostream>
using std::cout;

#include <boost/math/tools/roots.hpp>
using boost::math::tools::bisect;

#include "../source/risk_neutral.h"

#ifndef TOL
#define TOL 1e-3
#endif // !TOL

int main() {
    // *** Part 1 ***: Test with put option from market
    // Yahoo Finance (MSFT at 2024-02-16)
    const double target_growth = 0.12;  // Growth rate from Yahoo Finance (single stock)
    const double put_strike = 0.045; // Strike at chosen risk-free rate (U.S. Treasury Bill)

    // Corresponding put option rate from market corresponding to above
    const double put_price = 29.90 / 406.56;
    const double growth_variance = 0.2030; // Implied volatility from chosen option 

    cout << "Testing ROL calculations...with put option.\n";
    cout << "Put strike:            " << put_strike*100 << "%\n";
    cout << "Put price:             " << put_price*100 << "%\n";
    cout << "Target growth:         " << target_growth * 100 << "%\n";

    double i_trunc = get_returns_with_put(target_growth, growth_variance, put_price, put_strike);
    cout << "Growth with put:       " << i_trunc*100 << "%\n";
    double rol = get_min_capital_ROL(target_growth, put_price, i_trunc);
    cout << "ROL = " << rol*100 << "%\n\n";

    // *** Part 2 ***: Test with no put option
    cout << "Testing ROL calculations...without put option.\n";
    cout << "Using same target growth and risk-free rate.\n";
    double rol2 = get_min_capital_ROL(target_growth, 0.0, put_strike);
    cout << "ROL = " << rol2*100 << "%\n";

    cout << "Option ROL < Swap ROL? " << ((rol <= rol)?"Yes":"No") << "\n\n";

    // *** Part 3 ***: Implied target rate of return required on MSFT with given market ROL
    const double mkt_rol = 0.5 / 100;
    cout << "Using example market rate on line of " << mkt_rol*100 << "%, calculate require target rate of return on MSFT.\n";

    // Use bisect to find root
    std::pair<double, double> result = bisect(
        [&](double _tgt) { return get_min_capital_ROL(_tgt, put_price, i_trunc) - mkt_rol; },
        0.0,
        1.0,
        [](double l, double r) { return abs(l-r) < TOL; }
    );

    cout << "Implied target rate (before put):  " << (result.first + result.second)/2*100 << "%\n\n";
    return 0;
}
