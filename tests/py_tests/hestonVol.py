# Script to test fitting Heston model to live market data

# import numpy as np
# import pandas as pd

from src.pyvol import cpy_credit as cc
from src.pyvol.heston import get_call_information
from src.pyvol.heston import fit_heston as py_fit_Heston

# Get latest share information and collection of live traded options
print('Getting information on AAPL')
S0, vol_surf = get_call_information('AAPL')

print(f'Spot price of AAPL:     {S0:,.2f}\n')

# Remove any 0-options
vol_surf = vol_surf[vol_surf['price'] > 0]
vol_surf = vol_surf[vol_surf['maturity'] > 0]
print('Snippet of volatility surface data:')
print(vol_surf.head())

# Fit Heston model in C++
pyV = cc.fit_Heston(S0, vol_surf['strike'].values, vol_surf['rf'].values, vol_surf['maturity'].values, vol_surf['price'].values)

print('\nFit results (C++):')
print(f'Spot volatility (v0):           {pyV.v0:,.1%}')
print(f'Mean-rev rate (kappa):          {pyV.alpha:,.1f}')
print(f'Asymp mean var (theta):         {pyV.vTheta:,.1%}')
print(f'Vol of vol (sigma):             {pyV.vSig:,.1%}')
print(f'Markt price of vol (lambda):    {pyV.vLambda:,.1%}')
print(f'Price/vol corr (rho):           {pyV.rho:,.1%}')
print('\n')

# Fit Heston model in Python
print('Fit results (Python):')
resultHeston = py_fit_Heston(vol_surf, S0)
for key, item in resultHeston.items():
    print(f'{key}: {item}')
