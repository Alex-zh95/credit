# Script to test fitting Heston model to live market data

# import numpy as np
# import pandas as pd

from src.pyvol import cpy_credit as cc
from src.pyvol.heston import get_call_information

# Get latest share information and collection of live traded options
S0, vol_surf = get_call_information('AAPL')

# Remove any 0-options
vol_surf = vol_surf[vol_surf['price'] > 0]

# Fit Heston model
pyV = cc.fit_Hst(S0, vol_surf['strike'], vol_surf['maturity'], vol_surf['price'])

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
