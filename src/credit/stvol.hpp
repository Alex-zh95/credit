#ifndef STVOL_HPP
#define STVOL_HPP

#include <complex>
#include <memory>
#include <vector>

/* @Filename:        stvol.hpp
 * @Description:     Declares containers and objects for dealing with Heston model for volatility surface.
 */

// General Note: Move constructors implicitly generated as long as no copy constructor/assignment defined, other than delete, or user-generated destructor.

namespace StVol
{
    /* @Description:    Data structure to hold the parameters of underlying instrument.
     *
     * @Params:         double S0:          Spot/current price of underlying
     *                  double v0:          Spot/current volatility of underlying
     *                  double alpha:       Mean reversion rate
     *                  double vSig:        Volatility of volatility
     *                  double rho:         Correlation parameter (volatility to price)
     *                  double vTheta:      Long-run variance
     *                  double vLamdba:     Market price of volatility (risk premium)
     *                  double rf:          Risk-free interest rate
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

    /* @Description:    Model to encapsulate the Heston model.
     */
    class HestonCallMdl
    {
    public:
        // Main constructor
        HestonCallMdl(std::unique_ptr<Underlying> _underlying, double _K, double _t = 1.) : underlying(std::move(_underlying)), K(_K), t(_t), P(0.0)
        {}

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

        /* @Description: Main pricing function of the Heston model to calculate. */
        void calc_option_price();

    private:
        /* @Description: Define the characteristic function used in Heston valuation.
         *
         * @Params:     double phi:           Frequency-spectrum input to the characteristic function
         *
         * @Returns:    std::complex<double>: Result
         */
        std::complex<double> charFn(std::complex<double> phi);

        /* @Description: Defining the integrand needed for getting the price of the Heston model
         *
         * @Params:     double phi:           Frequency-spectrum input to the characteristic function
         *
         * @Returns:    std::complex<double>: Result
         */
        std::complex<double> integrand(double phi);

        std::unique_ptr<Underlying> underlying;
        double K;
        double t;
        double P;
    };

    /* @Description: Calibrate a Heston model to market data of traded call options. Method involves optimizing the square error function of expected and observed Heston call prices against market data.
     *
     * @Params:     std::vector<HestonCallMdl>& hMdls:  Vector of market-traded call options
     *              std::vector<double> market_prices:  Prices of those call options
     *              std::vector<double> initial_guess:  Initial estimate of values
     *
     * @Returns:    std::unique_ptr<Underlying>:        Optimized underlying parameters for the Heston model.
     */
    std::unique_ptr<Underlying> market_calibration(std::vector<HestonCallMdl>& hMdls, const std::vector<double> market_prices, std::vector<double> initial_guess);
} // namespace StVol

#endif // STVOL_HPP
