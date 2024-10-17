# Preamble for running script outside of REPLs: add parent directory to sys path
import sys
import os

script_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(script_dir))

import numpy as np
import yfinance as yf

from datetime import datetime
from src.pyvol.garch11 import GarchVol
from build import cpy_credit as cc


print('Testing garch11.py...')
now = datetime.today()
clip = datetime(now.year - 1, now.month, now.day)

# Assumptions
tick_name = 'ACGL'
TRADING_DAYS = 252
rf = 0.045  # Risk-free rate
infl = 0.03  # Inflation

print(f'Loading {tick_name} equity data...')
ticker = yf.Ticker(tick_name)
df = ticker.history(start=clip, end=now)
df.head()

# Use close values
mdl = GarchVol(X=df['Close'].values)

print('\nEmpirical daily return statistics:')
emp_sig, emp_mu = mdl.sigma, mdl.mu
print(f'Sigma:   {emp_sig:,.3f}%')
print(f'Mu:      {emp_mu:,.3f}%')

print('\nFitting GARCH...')
mdl.fit_garch()

fit_sig, fit_mu = np.mean(mdl.sigma), np.mean(mdl.mu)
print(f'Sigma:   {fit_sig:,.3f}%')
print(f'Mu:      {fit_mu:,.3f}%')

print('\nDiagnostics:')
test = mdl.diagnostics()
print(f'p-val:   {test[-1]}')

# Convert the stats to annual
vol = fit_sig * np.sqrt(TRADING_DAYS) / 100  # Use correct decimalization
drft = fit_mu * TRADING_DAYS / 100
print(f'\nAnnual volatility: {vol:,.3%}')
print(f'Annual drift:      {drft:,.3%}')

# For a basic FPT model, we analyze total revenue and total expense (we should in theory use NEP and Total Loss Adj. Expenses (Claims) but cannot programmatically obtain with yfinance)
income = ticker.financials.loc['Total Revenue'].iloc[0] / 1e9
outgo = ticker.financials.loc['Total Expenses'].iloc[0] / 1e9

print(f'\nTotal income (bns):    {income:.3f}')
print(f'Total outgo (bns):     {outgo:.3f}')
print(f'Combined ratio:        {outgo / income:,.3%}')

print('\nUsing the First-Passage Time model, calculate risk-neutral reinsurance default rate...')

asset_sig = cc.get_asset_volatility(income, vol, outgo, rf)
print(f'Obtain implied asset volatility: {asset_sig:,.3%}')

# Discounting to present used in probability calc
rn_rol = cc.get_fpt_default_probability(a0=income, rf=rf, sigma_a=asset_sig, L=outgo, gamma=infl) * np.exp(-rf)
print(f'Risk-neutral ROL:      {rn_rol:,.3%}')

# Risk-adjust via the Wang transform
rol = cc.wang_transform([rn_rol], sharpe_ratio=-drft)[0]
print(f'Risk-adj. ROL:         {rol:,.3%}')
