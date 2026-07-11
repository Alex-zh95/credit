#ifndef STVOL_HPP
#define STVOL_HPP

#include <complex>
#include <memory>
#include <vector>

#include "params.hpp"

// General Note: Move constructors implicitly generated as long as no copy constructor/assignment
// defined, other than delete, or user-generated destructor.

namespace StVol {
    /**
     * Description
     * -----------
     * Data structure to hold the parameters of an underlying instrument following
     * Heston (stochastic volatility) dynamics.
     *
     * Inherits the common contract/underlying parameters from StandardUnderlying
     * (S0, K, r, sigma, t, call, q), where:
     *      double S0:          Spot/current price of underlying
     *      double K:           Strike price
     *      double r:           Risk-free interest rate
     *      double t:           Duration of the contract (expiry)
     *
     * and extends with the stochastic volatility parameters:
     *
     * Params
     * ------
     * double v0:          Spot/current volatility of underlying
     * double alpha:       Mean reversion rate
     * double vSig:        Volatility of volatility
     * double rho:         Correlation parameter (volatility to price)
     * double vTheta:      Long-run variance
     * double vLamdba:     Market price of volatility (risk premium)
     */
    struct HestonUnderlying : StandardUnderlying {
        double v0;
        double alpha;
        double vSig;
        double rho;
        double vTheta;
        double vLambda;

        // Default constructor - all parameters zeroed (unit maturity, call from base)
        HestonUnderlying()
            : StandardUnderlying(), v0(0.), alpha(0.), vSig(0.), rho(0.), vTheta(0.), vLambda(0.) {}

        /**
         * Construct from base contract parameters plus a std::vector<double> of
         * volatility parameters, in the order: v0, alpha, vTheta, vSig, vLambda, rho
         */
        HestonUnderlying(const StandardUnderlying& params, const std::vector<double>& volParams)
            : StandardUnderlying(params), v0(volParams.at(0)), alpha(volParams.at(1)),
              vSig(volParams.at(3)), rho(volParams.at(5)), vTheta(volParams.at(2)),
              vLambda(volParams.at(4)) {}
    };

    class HestonCallMdl {
      public:
        // Main constructor - strike and maturity stored on the underlying (K, t)
        HestonCallMdl(std::unique_ptr<HestonUnderlying> _underlying, double _K, double _t = 1.)
            : underlying(std::move(_underlying)), P(0.0) {
            underlying->K = _K;
            underlying->t = _t;
        }

        // Alternative constructor - pass in a std::vector<double> for volatility parameters
        HestonCallMdl(const StandardUnderlying& params, const std::vector<double>& volParams)
            : underlying(std::make_unique<HestonUnderlying>(params, volParams)), P(0.0) {}

        // Copy constructor
        HestonCallMdl(const HestonCallMdl& other) : P(other.P) {
            if (other.underlying)
                underlying = std::make_unique<HestonUnderlying>(*other.underlying);
        }

        // Move constructor
        HestonCallMdl(HestonCallMdl&& other)
            : underlying(std::move(other.underlying)), P(other.P) {}

        // Copy assignment operator
        HestonCallMdl& operator=(const HestonCallMdl& other) {
            if (this != &other) {
                underlying = other.underlying
                                 ? std::make_unique<HestonUnderlying>(*other.underlying)
                                 : nullptr;
                P = other.P;
            }

            return *this;
        }

        // Move assignment operator
        HestonCallMdl& operator=(HestonCallMdl&& other) {
            if (this != &other) {
                underlying = std::move(other.underlying);
                P = other.P;
            }

            return *this;
        }

        // Getters and setters
        void set_strike(double _K) { underlying->K = _K; }
        void set_maturity(double _t) { underlying->t = _t; }
        void set_underlying(std::unique_ptr<HestonUnderlying> _underlying) {
            underlying = std::move(_underlying);
        }

        double get_strike() { return underlying->K; }
        double get_maturity() { return underlying->t; }
        double get_option_price() { return P; }

        // Read-only const access to the Underlying struct
        const HestonUnderlying& get_underlying() const { return *underlying; }

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

        std::unique_ptr<HestonUnderlying> underlying;
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
    std::unique_ptr<HestonUnderlying> fitHeston(double spot_price, std::vector<double> strikes,
                                                std::vector<double> r,
                                                std::vector<double> maturities,
                                                std::vector<double> market_prices,
                                                std::vector<double> trade_volumes);

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
     * double maturity:                Maturity of using a call option to model company financing
     * structure
     */

    std::unique_ptr<StVol::HestonUnderlying> HestonAssetVolatilityImplied(StVol::HestonCallMdl& mdl,
                                                                          double asset, double debt,
                                                                          double maturity);
} // namespace StVol

#endif // STVOL_HPP