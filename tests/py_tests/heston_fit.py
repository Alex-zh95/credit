# Script to test fitting Heston model to live market data

# import numpy as np
# import pandas as pd

from pyvol import cpy_credit as cc
from pyvol.market_data import get_call_information

# Get latest share information and collection of live traded options
symb = 'CCJ'
print(f'Getting information on {symb}')
S0, vol_surf = get_call_information(symb)

print(f'Spot price of {symb}:     {S0:,.2f}\n')

# Remove 0-prices
vol_surf = vol_surf[vol_surf['price'] > 0]
vol_surf = vol_surf.dropna()
print('Snippet of volatility surface data:')
print(vol_surf.head())

# Fit Heston model in C++
pyV = cc.fit_Heston(S0, vol_surf['strike'].values, vol_surf['rf'].values, vol_surf['maturity'].values, vol_surf['price'].values, vol_surf['volume'].values)

print('\nFit results (C++):')
print(f'Spot volatility (v0):           {pyV.v0:,.5f}')
print(f'Mean-rev rate (alpha):          {pyV.alpha:,.5f}')
print(f'Asymp mean var (theta):         {pyV.vTheta:,.5f}')
print(f'Vol of vol (sigma):             {pyV.vSig:,.5f}')
print(f'Markt price of vol (lambda):    {pyV.vLambda:,.5f}')
print(f'Price/vol corr (rho):           {pyV.rho:,.5f}')
print('\n')

# Suppose now we want to calculate the default probability of this firm with a
# given debt-to-asset ratio in 1 year's time
debt_ratio = 0.75
pDefault, pyAssetVol = cc.get_Heston_default_probability(pyV, 1.0, debt_ratio)
print(f'RN Prob of default:             {pDefault:,.5f}\n')

print('\nImplied asset volatility parameters (C++):')
print(f'Spot volatility (v0):           {pyAssetVol.v0:,.5f}')
print(f'Mean-rev rate (alpha):          {pyAssetVol.alpha:,.5f}')
print(f'Asymp mean var (theta):         {pyAssetVol.vTheta:,.5f}')
print(f'Vol of vol (sigma):             {pyAssetVol.vSig:,.5f}')
print(f'Markt price of vol (lambda):    {pyAssetVol.vLambda:,.5f}')
print(f'Price/vol corr (rho):           {pyAssetVol.rho:,.5f}')
print('\n')
