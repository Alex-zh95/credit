#ifndef RISK_NEUTRAL_H
#define RISK_NEUTRAL_H
/* @Filename:        risk_neutral.h
 * @Description:     Declares functions for dealing with risk-neutral probability space.
 */
#include <vector>

/* @Description:    Calculate the vanilla option price via the Black-Scholes equation, using the provided inputs:
 *
 * @Params:         double S0:          Current price of underlying
 *                  double K:           Strike price
 *                  double r:           Risk-free interest rate
 *                  double sigma:       Volatility of observed price
 *                  double t:           Duration of the option (expiry)
 *                  double q = 0:       Assumed constant and continuous dividend rate
 *                  bool call = true    Whether looking at call or put
 *
 * @Returns:        double& price       Price of option at time 0
 *                  double& Phi1        Delta of the option
 *                  double& Phi2        Risk-neutral probability that option is exercised
 */
void option_price(
    double S0,
    double K,
    double r,
    double sigma,
    double t,
    double& price,
    double& Phi1,
    double& Phi2,
    bool call = true,
    double q = 0
);

/* @Description:    Converts a given actual probability vector into risk-adjusted probability.
 *                  The Wang Transform is one such method of doing so, and uses the Sharpe Ratio to do so.
 *
 * @Params:         std::vector<double> P:      Input vector of probabilities
 *                  double sharpe_ratio:        Sharpe ratio
 *                  bool inverse = false:       Set to true to invert the transform
 *
 * @Returns:        std::vector<double> Q:      Output vector of probabilities
 */
std::vector<double> wang_transform(
    std::vector<double> P,
    double sharpe_ratio,
    bool inverse = false
);

/* @Description:    Derive asset volatility via equity volatility numerically.
 *                  Motivation for this is due to asset volatility not being observable.
 *
 * @Params:         std::vector<double> E:      Vector of stock values
 *                  double sigma_e:             Equity volatility
 *                  double L:                   Face value of debt
 *                  double r:                   Risk-free rate of return
 *                  double t:                   Debt duration, defaulted to 1
 *                  uint n_iter:                Maximum number of iterations, defaulted to 50
 *
 * @Returns:        double sigma_a:             Asset volatility
 */
double get_asset_volatility(
    std::vector<double> E,
    double sigma_e,
    double L,
    double r,
    double t = 1,
    unsigned int n_iter = 50
);

/* @Description:    Calculate risk-neutral probability of default.
 *
 * @Params:         double a0:              Current asset value
 *                  double mu_a:            Asset drift
 *                  double sigma_a:         Asset volatility
 *                  double L:               Debt face value
 *                  double t:               Duration of debt, defaulted to 1
 *
 * @Returns:        double p_default:       Probability of default
 */
double get_default_probability(
    double a0,
    double mu_a,
    double sigma_a,
    double L,
    double t = 1
);

#endif
