#ifndef PARAMS_HPP
#define PARAMS_HPP

/**
 * Description
 * -----------
 * Base data structure holding parameters common to underlying instruments and the
 * (vanilla) option contracts written on them. Specialized underlyings (e.g.
 * StVol::HestonUnderlying) inherit from this rather than redefining common symbols.
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
                       bool call_ = true, double q_ = 0.) noexcept
        : S0(S0_), K(K_), r(r_), sigma(sigma_), t(t_), call(call_), q(q_) {}

    StandardUnderlying() noexcept : StandardUnderlying(0., 0., 0., 0.) {}

    virtual ~StandardUnderlying() = default;
    StandardUnderlying(const StandardUnderlying&) = default;
    StandardUnderlying& operator=(const StandardUnderlying&) = default;
    StandardUnderlying(StandardUnderlying&&) = default;
    StandardUnderlying& operator=(StandardUnderlying&&) = default;
};

struct AssetDefaultParams {
    double a0;       // Current asset value
    double rf;       // Risk-free rate
    double sigma_a;  // Asset volatility
    double L;        // Debt face value
    double t;        // Maturity
    double q;        // Dividend rate

    AssetDefaultParams(double a0_, double rf_, double sigma_a_, double L_, double t_ = 1.,
                       double q_ = 0.) noexcept
        : a0(a0_), rf(rf_), sigma_a(sigma_a_), L(L_), t(t_), q(q_) {}

    AssetDefaultParams() noexcept : AssetDefaultParams(0., 0., 0., 0.) {}
};

#endif
