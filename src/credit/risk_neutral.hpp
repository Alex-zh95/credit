#ifndef RISK_NEUTRAL_H
#define RISK_NEUTRAL_H
/* @Filename:        risk_neutral.hpp
 * @Description:     Declares functions for dealing with risk-neutral probability space.
 */
#include <vector>
#include <tuple>

/* @Description:    Calculate the vanilla option price via the Black-Scholes equation, using the provided inputs:
 *
 * @Params:         double S0:          Current price of underlying
 *                  double K:           Strike price
 *                  double rf:          Risk-free interest rate
 *                  double sigma:       Volatility of observed price
 *                  double t = 1:       Duration of the option (expiry)
 *                  double q = 0:       Assumed constant and continuous dividend rate
 *                  bool call = true:   Whether looking at call or put, default is call
 *
 * @Returns:        std::tuple<double, double, double>, split into the following:
 *                      double price       Price of option at time 0
 *                      double Phi1        Delta of the option
 *                      double Phi2        Risk-neutral probability that option is exercised
 */
std::tuple<double, double, double> vanilla_option_price(
    const double S0,
    const double K,
    const double r,
    const double sigma,
    const double t = 1,
    const bool call = true,
    const double q = 0);

/* @Description:    Calculate the down-and-out call option price via the Black-Scholes equation, using the provided inputs:
 *
 * @Params:         double S0:          Current price of underlying
 *                  double K:           Strike price
 *                  double r:           Risk-free interest rate
 *                  double sigma:       Volatility of observed price
 *                  double t = 1:       Duration of the option (expiry)
 *                  double gamma = 0:   Growth rate of the strike price
 *                  double q = 0:       Dividend rate
 *
 * @Returns:        std::tuple<double, double>, split into the following:
 *
 *                      double c_do:    Price of down-and-out call at time 0
 *                      double c_di:    Price of down-and-in call at time 0
 */
std::tuple<double, double> fpt_call_price(
    const double S0,
    const double K,
    const double r,
    const double sigma,
    const double t = 1,
    const double gamma = 0,
    const double q = 0);

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
    const double sharpe_ratio,
    const bool inverse = false);

/* @Description:    Derive asset volatility via equity volatility numerically.
 *                  Motivation for this is due to asset volatility not being observable.
 *
 * @Params:         double E:                   Stock value
 *                  double sigma_e:             Equity volatility
 *                  double L:                   Face value of debt
 *                  double r:                   Risk-free rate of return
 *                  double t = 1:               Debt duration
 *                  const n_iter = 50:          Maximum number of iterations
 *
 * @Returns:        double sigma_a:             Asset volatility
 */
double get_asset_volatility(
    const double E,
    const double sigma_e,
    const double L,
    const double r,
    const double t = 1,
    const int n_iter = 50);

/* @Description:    Calculate risk-neutral probability of default under Merton model.
 *
 * @Params:         double a0:              Current asset value
 *                  double rf:              Asset drift
 *                  double sigma_a:         Asset volatility
 *                  double L:               Debt face value
 *                  double t = 1:           Duration of debt
 *
 * @Returns:        double p_default:       Probability of default
 */
double get_vanilla_default_probability(
    const double a0,
    const double rf,
    const double sigma_a,
    const double L,
    const double t = 1);

/* @Description:    Calculate risk-neutral probability of default under First Passage Time model
 *
 * @Params:         double a0:              Current asset value
 *                  double rf:              Asset drift
 *                  double sigma_a:         Asset volatility
 *                  double L:               Debt boundary value
 *                  double q = 0:           Dividend rate
 *                  double gamma = 1:       Debt growth rate
 *                  double t = 1:           Duration of debt
 *
 * @Returns:        double p_default:       Probability of default
 */
double get_fpt_default_probability(
    const double a0,
    const double rf,
    const double sigma_a,
    const double L,
    const double q = 0,
    const double gamma = 0,
    const double t = 1);

/* @Description:    Calculate minimum rate on line (ROL).
 *                  For the non-option method, set i=risk-free rate and p=0.
 *
 * @Params:         double y:               Target investment rate (yield),
 *                  double p:               Price of put option used to hedge against investment achieving less than risk-free return,
 *                  double i:               Expected investmet return with put protection.
 *
 * @Returns:        double min_ROL
 */
double get_min_capital_ROL(
    const double y,
    const double p,
    const double i);

/* @Description:        Calculate returns on investment with purchase of put option.
 *
 * @Params:             double y:           Original target investment rate (yield)
 *                      double y_var:       Investment volatility (e.g. implied from put price).
 *                      double put:         Market price of put option.
 *                      double r:           Risk-free rate.
 */
double get_returns_with_put(
    const double y,
    const double y_var,
    const double put,
    const double r);

/* @Description:    Calculate corresponding price of a put or call option given a call or put option via put-call parity.
 *
 * @Params:         double in_price:        Price of instrument to convert
 *                  double rf:              Risk-free interest rate
 *                  double S0:              Value of underlying instrument
 *                  double K:               Strike price
 *                  double t = 1:           Duration of contract
 *                  bool in_call = true:    Type of instrument, default to call
 * 
 * @Returns:        double out_price:       Price of converted instrument
 */
double put_call_parity(
    const double in_price,
    const double rf,
    const double S0,
    const double K,
    const double t = 1,
    bool in_call = true);

#endif // RISK_NEUTRAL_H