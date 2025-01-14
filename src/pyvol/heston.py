r'''
Demo file for fitting the Heston model for volatility surface.

The Heston model of equity is given by:

$$dS_t = rS_t dt + \sqrt{v_t} S_t dW_t$$
$$dv_t = \kappa (\tehta - v_t) dt + \sigma \sqrt{v_t} dB_t$$
$$\rho dt = dW_t dB_t$$

under risk-neutral valuation. Here,

- $S_t$ is the equity spot price,
- $v_t$ is the spot variance,
- $K$ is strike price,
- $W_t, B_t$ are standard Brownian motions,
- $r_t$ is the risk-free interest rate,
- $\kappa$ is the mean reversion rate,
- $\theta$ is the long-run variance,
- $\sigma$ is the volatility of variance,
- $\rho$ is the correlation parameter.

To convert from real-world to risk-neutral, we require $\lambda$, the variance risk premium:

$$\kappa = \kappa' + \lambda$$
$$\theta = \frac{\kappa' \theta'}{\kappa' + \lambda}$$

where $\kappa', \theta'$ correspond to mean reversion rate and long-run variance but under real-world measure.

Heston's model is solved through the use of characteristic functions, $\varphi_1, \varphi_2$. The solution itself for a vanilla European call option is given by:

$$C(t) = \frac{1}{2}(S_0 - Ke^{-rt}) + \frac{1}{\pi} \int_0^{\inf} \Re\left\[ e^{rt} \frac{\varphi (\theta - i)}{i \phi K^{i\phi}} - K \frac{\varphi(\theta)}{i\phi K^{i\theta}} \right\] d\phi$$

where $\varphi$ is a simplification for $\varphi_1, \varphi_2$ through some change of variables.
'''

import numpy as np
import pandas as pd
from scipy.integrate import quad
from scipy.optimize import minimize
from datetime import datetime
import yfinance as yf

from src.pyvol.yields import get_latest_yields


def heston_char(phi, S0, v0, kappa, theta, sigma, rho, lambd, tau, r):
    '''
    Define the Heston characteristic function using the parameters as provided.
    '''
    # Define some constants to simplify
    a = kappa * theta
    b = kappa + lambd
    rspi = rho * sigma * phi * 1j
    d = np.sqrt((rspi - b)**2 + (phi * 1j + phi**2) * sigma**2)
    g = (b - rspi + d) / (b - rspi - d)

    exp1 = np.exp(r * phi * 1j * tau)
    term2 = S0**(phi * 1j) * ((1 - g * np.exp(d * tau)) / (1 - g))**(-2 * a / sigma**2)
    exp2 = np.exp(a * tau * (b - rspi + d) / sigma**2 + v0 * (b - rspi + d) * ((1 - np.exp(d * tau)) / (1 - g * np.exp(d * tau))) / sigma**2)

    return exp1 * term2 * exp2


def heston_integrand(phi, S0, v0, kappa, theta, sigma, rho, lambd, tau, r, K):
    char1_params = (phi - 1j, S0, v0, kappa, theta, sigma, rho, lambd, tau, r)
    char2_params = (phi, S0, v0, kappa, theta, sigma, rho, lambd, tau, r)

    numerator = np.exp(r * tau) * heston_char(*char1_params) - K * heston_char(*char2_params)
    denominator = 1j * phi * K**(1j * phi)

    return numerator / denominator


def heston_price(S0, K, v0, kappa, theta, sigma, rho, lambd, tau, r):
    args = (S0, v0, kappa, theta, sigma, rho, lambd, tau, r, K)

    integrated, _ = np.real(quad(heston_integrand, 0, 100, args=args))
    return (S0 - K * np.exp(-r * tau)) / 2 + integrated / np.pi


def get_call_information(ticker_symb: str) -> tuple[float, pd.DataFrame]:
    '''
    Step 0: Retrieve options from Yahoo finance into a dataframe, containing:
    - expiration date,
    - strike,
    - price = (bid+ask)/2

    Step 1: Generate a volatility surface dataframe with index=maturities, cols=strikes
    Step 2: Melt the above to get a long list of options with cols maturities, strikes, price
    Step 3: Generate a risk-free rate for each maturity via a yield curve (extra col)
    '''

    # Step 0: Handle to ticker object
    tk = yf.Ticker(ticker_symb)

    # Get most recent price
    S0 = tk.history()['Close'].iloc[-1]

    # Get option data
    exps = tk.options  # Read possible expiration dates
    options = []
    for e in exps:
        opt = tk.option_chain(e)
        opt = opt.calls
        opt['expiry_date'] = datetime.strptime(e, '%Y-%m-%d')
        options.append(opt)

    options = pd.concat(options, axis=0).reset_index()

    options['maturity'] = (options['expiry_date'] - datetime.today()).dt.days / 365.25

    options[['bid', 'ask', 'strike']] = options[['bid', 'ask', 'strike']].apply(pd.to_numeric)

    options['price'] = options[['bid', 'ask']].apply(np.mean, axis=1)

    # Steps 1 and 2: Volatility surface
    vol_surface = options[['maturity', 'strike', 'price']]

    # Step 3: Apply a yield curve to generate the risk-free rate applicable for each maturity
    y_curve = get_latest_yields()
    vol_surface['rf'] = vol_surface['maturity'].apply(y_curve)

    return S0, vol_surface


def fit_heston(vol_surface: pd.DataFrame, S0: float, rate_col_str='rf', strike_col_str='strike', maturity_col_str='maturity', price_col_str='price') -> dict:
    r = vol_surface[rate_col_str].values
    K = vol_surface[strike_col_str].values
    tau = vol_surface[maturity_col_str].values
    P = vol_surface[price_col_str].values

    # Set up initial values and lower and upper bounds for each parameter
    params = {
        'v0': {'x0': 0.1, 'lbub': [1e-3, 0.1]},
        'kappa': {'x0': 3, 'lbub': [1e-3, 5]},
        'theta': {'x0': 0.05, 'lbub': [1e-3, 0.1]},
        'sigma': {'x0': 0.3, 'lbub': [1e-2, 1]},
        'rho': {'x0': -0.8, 'lbub': [-1, 0]},
        'lambd': {'x0': 0.03, 'lbub': [-1, 1]}
    }

    x0 = [param['x0'] for _, param in params.items()]
    bnds = [param['lbub'] for _, param in params.items()]

    def _square_error(x):
        v0, kappa, theta, sigma, rho, lambd = [param for param in x]

        err = 0.0

        for i in range(vol_surface.shape[0]):
            err += (P[i] - heston_price(S0, K[i], v0, kappa, theta, sigma, rho, lambd, tau[i], r[i]))**2

        return err / vol_surface.shape[0]

    result = minimize(_square_error, x0, tol=1e-3, method='SLSQP', options={'maxiter': 1e4}, bounds=bnds)

    v0, kappa, theta, sigma, rho, lambd = [param for param in result.x]

    out = {
        'v0': v0,
        'kappa': kappa,
        'theta': theta,
        'sigma': sigma,
        'rho': rho,
        'lambd': lambd
    }

    return out


# %% TEST FUNCTIONS
def pricing_test() -> int:
    '''
    Test to see if implementation yields expected price. Returns 0 for no errors.
    '''
    S0 = 100.
    K = 100.
    v0 = 0.1
    r = 0.03
    kappa = 1.5768
    theta = 0.0398
    sigma = 0.3
    lambd = 0.575
    rho = -0.5711
    tau = 1.
    expected = 11.5403

    C = heston_price(S0, K, v0, kappa, theta, sigma, rho, lambd, tau, r)

    return 0 if np.abs(C - expected) / expected < 1e-3 else 1


def heston_test():
    '''
    Tests to see run-through of the Heston fitting algorithm on new data.
    '''
    S0, vol_surf = get_call_information('AAPL')
    vol_surf = vol_surf[vol_surf['price'] > 0]
    resultHeston = fit_heston(vol_surf, S0)

    for key, item in resultHeston.items():
        print(f'{key}: {item}')


if __name__ == '__main__':
    print('Testing Heston pricer...')
    err = pricing_test()

    print(f'Error code:     {err}\n\n')

    print('Testing Heston fitting via Python...')
    heston_test()
