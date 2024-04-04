# %% [markdown]
# # Merton-based reinsurance pricing example
#
# In this example, we look at deriving minimum rate on line (ROL) for AIG. Unlike regular insurance pricing, which looks at loss cost, we will be focusing effectively on capacity only (i.e. very high layer excesses, where loss costs cannot be modeled). We consider details from AIG's financials and investor activities.
#
# We assume that the company's insurance funds is comprised solely by capital and reserves. In equation form,
#
# $$F = C + R$$
#
# In one year's time, the result of the policy to the company is $\max(F-R, 0)$. In other words, we can model the value of the funds as a call on the company's funds with strike price set to the value of the reserves (the Merton model).
#
# Assumptions: $F$ can be modeled through a Geometric Brownian Motion and that we have certainty on reserves at the start of the year.
#
# The rate on line would effectively be the probability of default, calculated from the Merton model.
#
# ## Balance Sheet Data
#
# - Page 65 of the AIG's 2022 Annual Statement contains information on the return on equity.
# - Page 268 of Annual Statement contains insurance information. We take the data for "General Insurance".

# %%
import pandas as pd
import yfinance as yf
import datetime as dt
import numpy as np

from build import cpy_credit as cc

notebook_creation_date = dt.datetime.today()

financials = {
    'Net Premiums': 25.512,
    'Net Investment Income': 2.382,
    'Losses': 15.407,
    'Other Expenses': 3.533 + 4.352,  # Operating expenses
    'ROE': 0.21,
    'Symbol': 'AIG'
}

# %% [markdown]
#
# In order to evaluate for a "real world" probability of default, we use the the return on equity as the drift parameter. For a risk-neutral view, we would just use the risk-free rate, which here is the U.S. Treasury Bill rate.

# %%
drift = financials['ROE']
print(f'Growth parameter assumed:               {drift:,.3%}')

rf = 4.75 / 100  # 52-week U.S. T-bill discount rate
print(f'U.S. T-bill (risk-free) rate:           {rf:,.3%}')

# Projected interest payout
capital = financials['Net Premiums'] * (1 + drift) + financials['Net Investment Income'] * (1 + rf)
reserve = (financials['Losses'] + financials['Other Expenses']) * (1 + drift)

print(f'Total capital:                          {capital:,.3f}')
print(f'Total reserve:                          {reserve:,.3f}')

# %% [markdown]
# ## Volatility
#
# We cannot observe asset or fund volatilities so instead, we will imply it from equity information. We will also download option data with expiry closest to end-of-year and closest to at-the-money.

# %%
ticker = yf.Ticker(financials['Symbol'])

# Vector of daily close prices
equity = ticker.history(
    start=dt.datetime(notebook_creation_date.year - 1, 1, 1),
    end=notebook_creation_date,
)['Close'].values

print(f'Current price (close):                  {equity[-1]:,.2f}')


def select_equity_volatility(ticker, equity_price_today, end_date) -> float:
    expiry_dates = ticker.options

    options = pd.DataFrame()

    # Obtain all the put options
    for T in expiry_dates:
        cur_T_options = ticker.option_chain(T)
        options = pd.concat([options, cur_T_options.puts], ignore_index=True)
        options['expiry'] = T

    # Odd issue that yields wrong expiration dates so add 1 day to correct
    options['expiry'] = pd.to_datetime(options['expiry']) + dt.timedelta(days=1)
    options['duration'] = np.abs(((options['expiry']) - end_date).dt.days / 365)

    # Also look for only those options in the money
    options['strike'] = options['strike'].apply(pd.to_numeric)
    options = options[options['strike'] >= equity_price_today]

    # After sorting, return the average of the top 5 options
    options = options.sort_values(by=['duration', 'strike'], ascending=True)
    return np.mean(options['impliedVolatility'].iloc[:5])


# Look for the option closest to this price with expiry arround 1 year's time, taking average volatilities
equity_volatility = select_equity_volatility(ticker, equity[-1], dt.datetime(
    notebook_creation_date.year + 1,
    notebook_creation_date.month,
    notebook_creation_date.day,
))

print(f'Implied equity volatility from option:  {equity_volatility:,.3%}')

# Attempt to calculate asset volatility if it fails, fall back to equity volatility.
try:
    asset_volatility = cc.get_vanilla_asset_volatility(equity, equity_volatility, reserve, rf)
    print(f'Implied asset volatility:               {asset_volatility:,.3%}')
except RuntimeError:
    print("No convergence to asset volatility. Using equity volatility as proxy...")
    asset_volatility = equity_volatility

# %% [markdown]
# ## Rate on Line Result
#
# We calculate default probability, $q$, using both to see sensitivity in choice of volatility. We could evaluate the reinsurance product as having this default probability of paying out a large limit or non-default probability of paying out 0.
#
# Hence the fair value of a \$ 1 limit excess of loss is the expected discounted value of the payout:
#
# $$V = e^{-rT} E(V_T|\mathcal{F}_0) = e^{-r} q$$
#
# where $T$ is the average duration of the contract and $\mathcal{F}_t$ is the filtration at time $t$. Taking this to be 1 for most insurance classes, we have:

# %%
q = cc.get_vanilla_default_probability(
    capital,
    drift,
    asset_volatility,
    reserve
)

print(f'Default probability:                    {q:,.3%}')

rol = np.exp(-rf) * q
print(f'Implied rate on line:                   {rol:,.3%}')
