import pandas as pd
import numpy as np
from src.pyvol.garch11 import GarchVol


# %%
print('Testing garch11.py...')

print('Loading file "eq.csv"...')
df = pd.read_csv('tests/py_tests/eq.csv')
df.head()

# Use close values
mdl = GarchVol(X=df['Close'].values)

print('\nEmpirical daily return statistics:')
emp_sig, emp_mu = mdl.sigma, mdl.mu
print(f'Sigma:   {emp_sig}')
print(f'Mu:      {emp_mu}')

print('\nFitting GARCH...')
mdl.fit_garch()

fit_sig, fit_mu = np.mean(mdl.sigma), np.mean(mdl.mu)
print(f'Sigma:   {fit_sig}')
print(f'Mu:      {fit_mu}')

print('\nDiagnostics:')
test = mdl.diagnostics()
print(f'p-val:   {test[-1]}')

# Convert the stats to annual
TRADING_DAYS = 252
vol = fit_sig * np.sqrt(TRADING_DAYS)
print(f'\nAnnual volatility: {vol}')
