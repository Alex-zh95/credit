'''
Collect the most recent yield curves with smoothing.
'''

from datetime import datetime, timedelta
import pandas as pd
import numpy as np
import requests
import io
from bs4 import BeautifulSoup

from nelson_siegel_svensson.calibrate import calibrate_ns_ols
import yfinance as yf


def _extract_maturity_from_string(maturity_string: str) -> float | str:
    '''
    Internal fn to convert strings like "1 Mo", "2 Mo" to 1/12, 2/12 etc. and "1 Yr", "5 Yr" to 1., 5. etc.

    On failure, return the string
    '''
    t = maturity_string.split(' ')[0]

    if 'Mo' in maturity_string:
        t = float(t) / 12
        return t
    elif 'Yr' in maturity_string:
        t = float(t)
        return t
    else:
        return maturity_string


def get_latest_yields():
    '''
    Attain from the US Treasury the latest month's yields and return a smoothing.

    Returns
    -------
    NelsonSiegelCurve: Takes input maturity (in years) to return the yield.
    '''
    url = 'https://home.treasury.gov/resource-center/data-chart-center/interest-rates/TextView?type=daily_treasury_yield_curve&field_tdr_date_value_month='

    # Find the latest date
    now = datetime.now()
    latest = datetime.strftime(now, '%Y%m')

    target_url = f'{url}{latest}'

    # Using this target URL, download the data
    response = requests.get(target_url)

    if response.status_code != 200:
        raise Exception('Access problem')

    soup = BeautifulSoup(response.text, 'html.parser')

    table = soup.find('table')
    if table is None:
        raise Exception("Cannot find table on site")

    df = pd.read_html(io.StringIO(str(table)))[0]

    # Cols to keep
    maturity_cols = df.columns
    maturity_cols = maturity_cols[(maturity_cols.str.contains('Yr')) | maturity_cols.str.contains('Mo')]

    df = df[maturity_cols]

    renamer = {maturity_cols[i]: _extract_maturity_from_string(maturity_cols[i]) for i in range(len(maturity_cols))}
    df = df.rename(columns=renamer)

    # To generate a smoothed yield curve, take the average over each of the maturities and then regress
    Y = df.mean(axis=0) / 100

    curve, sts = calibrate_ns_ols(np.array(Y.index, dtype=float), Y.values, tau0=1.0)

    if sts.status > 0:
        Warning(f'Message: {sts.message}')

    # Curve object can be used to provide smoothness
    return curve


def get_call_information(ticker_symb: str) -> tuple[float, pd.DataFrame]:
    '''
    Step 0: Retrieve options from Yahoo finance into a dataframe, containing:
    - expiration date,
    - strike,
    - price = (bid+ask)/2

    Step 1: Generate a volatility surface dataframe with index=maturities, cols=strikes
    Step 2: Melt the above to get a long list of options with cols maturities, strikes, price
    Step 3: Generate a risk-free rate for each maturity via a yield curve (extra col)

    Params
    ------
    ticker_symb: str
        Selected ticker string for Yahoo Finance.

    Returns
    -------
    S0: float
        Spot value of ticker.

    vol_surface: pd.DataFrame:
        Options data with columns ('maturity', 'strike', 'price', 'rf', 'volume').
    '''

    # Step 0: Handle to ticker object
    tk = yf.Ticker(ticker_symb)

    # Get most recent price - limit the data download to 5 days
    S0 = tk.history(period='5d')['Close'].iloc[-1]

    # Steps 1 and 2: Get option data and build volatility surface data
    exps = tk.options  # Read possible expiration dates
    vol_surface = []
    for e in exps:
        opt = tk.option_chain(e)
        opt = opt.calls

        # Expiry dates appear to be offset by 1 day - noticed in testing and also see
        # https://medium.com/@txlian13/webscrapping-options-data-with-python-and-yfinance-e4deb0124613
        opt['expiry_date'] = datetime.strptime(e, '%Y-%m-%d') + timedelta(days=1)
        vol_surface.append(opt)

    vol_surface = pd.concat(vol_surface, axis=0).reset_index()

    vol_surface['maturity'] = (vol_surface['expiry_date'] - datetime.today()).dt.days / 365.25
    vol_surface[['bid', 'ask', 'strike']] = vol_surface[['bid', 'ask', 'strike']].apply(pd.to_numeric)
    vol_surface['price'] = vol_surface[['bid', 'ask']].apply(np.mean, axis=1)

    # Step 3: Apply a yield curve to generate the risk-free rate applicable for each maturity
    y_curve = get_latest_yields()
    vol_surface['rf'] = vol_surface['maturity'].apply(y_curve)

    return S0, vol_surface[['maturity', 'strike', 'price', 'rf', 'volume', 'impliedVolatility', 'inTheMoney']]
