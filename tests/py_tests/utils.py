"""
Utility file to store Python runtime functions.
"""

import yfinance as yf
import pandas as pd
import numpy as np
import datetime as dt


def select_equity_volatility(ticker, equity_price_today, end_date) -> float:
    """
    Download equity volatility based on the yfinance ticker object given

    @Params:    ticker:             yfinance.ticker     Initialized ticker object
                equity_price_today: float               Current price of equity
                end_date:           dt.datetime         Target expiry date

    @Returns:   Mean implied volatility from options with selected expiry closest to the money.
    """
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
