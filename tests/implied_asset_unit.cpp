#include <iostream>
#include <vector>
#include "../source/risk_neutral.hpp"
#include <cmath>

std::vector<double> m_equity = {
    71.050003,
    71.730003,
    69.230003,
    59.07,
    58.529999,
    55.43,
    57.439999,
    60.009998,
    61.810001,
    58.919998,
    58.709999,
    51.82,
    51.700001,
    50.25,
    50.380001,
    53.759998,
    55.16,
    53.330002,
    56.529999,
    57.41,
    59.0,
    61.98,
    62.970001,
    64.18,
    63.57,
    61.950001,
    62.529999,
    63.77,
    63.990002,
    66.360001,
    63.060001,
    60.259998,
    62.18,
    63.939999,
    59.209999,
    60.009998,
    60.049999,
    62.41,
    62.25,
    63.16,
    63.540001,
    66.230003,
    65.220001,
    66.290001,
    68.919998,
    68.830002,
    69.5,
    70.360001,
    65.809998,
    68.550003,
};

const double rf = 4.750/100;
const double implied_equity_vol = 21.790/100;
const double reserve = 50.612;

const double expected_asset_vol = 12.787/100;

int main() {
    std::cout << "Testing get_asset_volatility...\n";
    std::cout << "Expected asset volatility:        " << expected_asset_vol << "\n";
    double asset_vol = get_vanilla_asset_volatility(m_equity, implied_equity_vol, reserve, rf);
    std::cout << "Actual asset volatility:          " << asset_vol << "\n";

    std::cout << "Test result?                      " << 
        ((std::abs(asset_vol - expected_asset_vol) < 1e-3)?"Pass":"Fail") << 
        std::endl;

    return 0;
}
