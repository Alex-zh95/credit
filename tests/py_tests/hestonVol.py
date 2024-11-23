# Script to test the Heston models

import numpy as np
from src.pyvol import cpy_credit as cc

from multiprocessing import Process, Queue

use_multiprocessing: bool = False

# Put a wrapper around the C++ function for fit_Hst for multiprocessing purposes
def py_fit_Hst(underlying, strikes, optn_prices, maturities, result_queue):
    result = cc.fit_Hst(underlying, strikes, maturities, optn_prices)
    result_queue.put(result)


# Construct some options data (AAPL as ref)
spot_price = 229.87

options_data = {
    'rf': [0.030, 0.035, 0.032, 0.033, 0.031, 0.029],
    'prices': [34.25, 31.20, 28.15, 24.40, 21.10],
    'strikes': [215., 220., 225., 230., 235.],
    'iv': [0.3256, 0.3197, 0.3121, 0.2945, 0.2806]
}

print('Testing fitting a Heston model and non-exercise probabilities...')

print('Using the following fictional market data:')
print(f'Spot price:     {spot_price:,.1f}')
print('Options data in dict form:')
print(options_data)
print('\n')

# Fit a Heston model to the above data - need to construct an Underlying object with initial guesses
pyU = cc.Underlying()
pyU.S0 = spot_price
pyU.v0 = 0.1
pyU.alpha = 1.5
pyU.vSig = 0.3           
pyU.rho = 0.05          # Assume very weak positive corr
pyU.vTheta = 0.1
pyU.rf = 0.03
pyU.vLambda = 0.3118    # 52-week change in market price

# pyU.v0 = 0.1
# pyU.alpha = 1.5
# pyU.vSig = 0.3           
# pyU.rho = 0.05          # Assume very weak positive corr
# pyU.vTheta = 0.1
# pyU.rf = 0.03
# pyU.vLambda = 0.3

# A fitted Heston yields another Underlying with the parameters
print('Fitting...')

if use_multiprocessing:
    result_queue = Queue()
    proc = Process(target=py_fit_Hst, args=(pyU, options_data['strikes'], np.repeat(1., 5), options_data['prices'], result_queue))
    proc.start()
    proc.join(timeout=5.0)

    if proc.is_alive():
        proc.terminate()
        proc.join()
        raise TimeoutError("Fitting Heston model time out. Terminating...")

    pyV = result_queue.get()
else:
    pyV = cc.fit_Hst(pyU, options_data['strikes'], np.repeat(1., 5), options_data['prices'])

print('\nFit results:')
print(f'Spot price:         {pyV.S0:,.1f}')
print(f'Spot volatility:    {pyV.v0:,.1%}')
print(f'Risk-free rate:     {pyV.rf:,.1%}')
print(f'Mean-rev rate:      {pyV.alpha:,.1f}')
print(f'Asymp mean var:     {pyV.vTheta:,.1%}')
print(f'Vol of vol:         {pyV.vSig:,.1%}')
print(f'Markt price of vol: {pyV.vLambda:,.1%}')
print(f'Price/vol corr:     {pyV.rho:,.1%}')
print('\n')

# Calculate the probability of not exercising the call option at the following strikes
strikes = spot_price * np.array([0.8, 0.9, 0.95, 0.99, 1.1])

for s in strikes:
    pOut = cc.get_Hst_default_probability(pyV, s)
    print(f'Probability out-of-money for spot {s:,.1f}: {pOut:.3%}')
