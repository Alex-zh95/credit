#ifndef STVOL_HPP
#define STVOL_HPP

#include <complex>
#include <memory>
#include <vector>

// General Note: Move constructors implicitly generated as long as no copy constructor/assignment defined, other than delete, or user-generated destructor.

namespace StVol
{
/** 
 * Description
 * -----------
 * Data structure to hold the parameters of underlying instrument.
 *
 * Params
 * ------
 * double S0:          Spot/current price of underlying
 * double v0:          Spot/current volatility of underlying
 * double alpha:       Mean reversion rate
 * double vSig:        Volatility of volatility
 * double rho:         Correlation parameter (volatility to price)
 * double vTheta:      Long-run variance
 * double vLamdba:     Market price of volatility (risk premium)
 * double rf:          Risk-free interest rate
 */
struct Underlying
{
    double S0;
    double v0;
    double alpha;
    double vSig;
    double rho;
    double vTheta;
    double vLambda;
    double rf;
};


class HestonCallMdl
{
public:
    // Main constructor
    HestonCallMdl(std::unique_ptr<Underlying> _underlying, double _K, double _t = 1.) : underlying(std::move(_underlying)), K(_K), t(_t), P(0.0)
    {}

    // Alternative constructor - pass in a std::vector<double> for volatility parameters
    HestonCallMdl(double _S0, const std::vector<double> &volParams, double _rf, double _K, double _t = 1.): K(_K), t(_t), P(0.0)
    {
        underlying = std::make_unique<Underlying>();
        underlying->S0 = _S0;
        underlying->v0 = volParams.at(0);
        underlying->alpha = volParams.at(1);
        underlying->vTheta = volParams.at(2);
        underlying->vSig = volParams.at(3);
        underlying->vLambda = volParams.at(4);
        underlying->rho = volParams.at(5);
        underlying->rf = _rf;
    }

    // Copy constructor
    HestonCallMdl(const HestonCallMdl& other)
    {
        if (other.underlying)
            underlying = std::make_unique<Underlying>(*other.underlying);
        t = other.t;
        P = other.P;
        K = other.K;
    }

    // Move constructor
    HestonCallMdl(HestonCallMdl&& other)
    {
        if (this != &other)
        {
            underlying = std::move(other.underlying);
            t = other.t;
            P = other.P;
            K = other.K;
        }
    }

    // Copy assignment operator
    HestonCallMdl& operator=(const HestonCallMdl& other)
    {
        if (this != &other)
        {
            underlying = std::make_unique<Underlying>(*other.underlying);
            t = other.t;
            P = other.P;
            K = other.K;
        }

        return *this;
    }

    // Move assignment operator
    HestonCallMdl& operator=(HestonCallMdl&& other)
    {
        if (this != &other)
        {
            underlying = std::move(other.underlying);
            t = other.t;
            P = other.P;
            K = other.K;
        }

        return *this;
    }


    // Getters and setters
    void set_strike(double _K) { K = _K; }
    void set_maturity(double _t) { t = _t; }
    void set_underlying(std::unique_ptr<Underlying> _underlying) { underlying = std::move(_underlying); }

    double get_strike() { return K; }
    double get_maturity() { return t; }
    double get_option_price() { return P; }

    // Read-only const access to the Underlying struct
    const Underlying& get_underlying() const { return *underlying; }

    // Get risk-neutral probability of exercise (asset > strike)
    double get_rn_exercise_probability();

    // Get delta of European call option.
    double get_delta();

    // Description: Main pricing function of the Heston model to calculate.
    void calc_option_price();

private:
    /**
     * Description
     * -----------
     * Define the characteristic function used in Heston valuation.
     *
     * Params
     * ------
     * double phi:           Frequency-spectrum input to the characteristic function
     *
     * Returns
     * -------
     * std::complex<double>: Result
     */
    std::complex<double> charFn(std::complex<double> phi);

    /** 
     * Description
     * -----------
     * Defining the integrand needed for getting the price of the Heston model
     *
     * Params
     * ------
     * double phi:           Frequency-spectrum input to the characteristic function
     *
     * Returns
     * -------
     * std::complex<double>: Result
     */
    std::complex<double> integrand(double phi);

    std::unique_ptr<Underlying> underlying;
    double K;
    double t;
    double P;
};

/** 
 * Description
 * -----------
 * Calibrate a Heston model to market data of traded call options by minimizing the
 * square error between Heston-implied prices against market prices.
 *  
 * Params
 * ------
 * double spot_price:                  Current price of asset
 * std::vector<double> strikes:        Vector of strike prices for options
 * std::vector<double> r:              Vector of risk-free rates
 * std::vector<double> maturities:     Vector of option durations
 * std::vector<double> market_prices:  Prices of those call options
 * std::vector<double> trade_volumes:  Volume of options traded.
 *
 * Returns:    std::unique_ptr<Underlying>:        Optimized underlying parameters
 */
std::unique_ptr<Underlying> fitHeston(
    double spot_price,
    std::vector<double> strikes,
    std::vector<double> r,
    std::vector<double> maturities,
    std::vector<double> market_prices,
    std::vector<double> trade_volumes
);

/**
 * Description
 * -----------
 * Implied asset volatility parameters, using idea of Merton structural model.
 * 
 * Params
 * ------
 * HestonCallMdl& mdl:             Reference to Heston model object
 * double asset:                   Asset value of underlying company
 * double debt:                    Debt value of the underlying company
 * double maturity:                Maturity of using a call option to model company financing structure
 */

std::unique_ptr<StVol::Underlying> HestonAssetVolatilityImplied(
    StVol::HestonCallMdl& mdl,
    double asset,
    double debt,
    double maturity
);
} // namespace StVol

#endif // STVOL_HPP
