#ifndef PARAMS_HPP
#define PARAMS_HPP

/**
 * Description
 * -----------
 * Base data structure holding parameters common to underlying instruments and the
 * (vanilla) option contracts written on them. Specialized underlyings (e.g.
 * StVol::HestonUnderlying, AssetDefaultParams) inherit from this rather than
 * redefining common symbols.
 *
 * Params
 * ------
 * double S0:          Spot/current price of underlying
 * double K:           Strike price
 * double r:           Risk-free interest rate
 * double sigma:       Volatility of observed price
 * double t:           Duration of the contract (expiry)
 * bool call:          Whether the contract is a call or a put
 * double q:           Assumed constant and continuous dividend rate
 */
struct StandardUnderlying {
    double S0;
    double K;
    double r;
    double sigma;
    double t;
    bool call;
    double q;

    StandardUnderlying(double S0_, double K_, double r_, double sigma_, double t_ = 1.,
                       bool call_ = true, double q_ = 0.)
        : S0(S0_), K(K_), r(r_), sigma(sigma_), t(t_), call(call_), q(q_) {}

    // Default constructor: prices, rates and volatility zeroed, unit maturity, call contract
    StandardUnderlying() : StandardUnderlying(0., 0., 0., 0.) {}
};

struct AssetDefaultParams : StandardUnderlying {
    // a0 maps to S0, rf maps to r, sigma_a maps to sigma, L maps to K
    AssetDefaultParams(double a0_, double rf_, double sigma_a_, double L_, double t_ = 1.,
                       double q_ = 0.)
        : StandardUnderlying(a0_, L_, rf_, sigma_a_, t_, true, q_) {}
};

#endif
