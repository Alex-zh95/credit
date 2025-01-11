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
from scipy.integrate import quad
from scipy.optimize import minimize


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
    exp2 = np.exp(a * tau * (b - rspi + d) / sigma**2) + v0 * (b - rspi + d) * ((1 - np.exp(d * tau)) / (1 - g * np.exp(d * tau))) / sigma**2

    return exp1 * term2 * exp2


def heston_integrand(phi, S0, v0, kappa, theta, sigma, rho, lambd, tau, r, K):
    char1_params = (phi - 1j, S0, v0, kappa, theta, sigma, rho, lambd, tau, r)
    char2_params = (phi, S0, v0, kappa, theta, sigma, rho, lambd, tau, r)

    numerator = np.exp(r * tau) * heston_char(*char1_params) - K * heston_char(*char2_params)
    denominator = 1j * phi * K**(1j * phi)

    return numerator / denominator


def heston_price(S0, K, v0, kappa, theta, sigma, rho, lambd, tau, r):
    args = (S0, v0, kappa, theta, sigma, rho, lambd, tau, r, K)

    integrated = np.real(quad(heston_integrand, 0, 100, args=args))
    return (S0 - K * np.exp(-r * tau)) / 2 + integrated / np.pi


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


if __name__ == '__main__':
    print('Testing Heston pricer...')
    err = pricing_test()

    print(f'Error code:     {err}')
