'''
Collect the most recent yield curves with smoothing.
'''

from datetime import datetime
import pandas as pd
import numpy as np
import requests
from bs4 import BeautifulSoup

from nelson_siegel_svensson.calibrate import calibrate_ns_ols


def maturity_string_to_years(maturity_string: str) -> float | str:
    '''
    Convert strings like "1 Mo", "2 Mo" to 1/12, 2/12 etc. and "1 Yr", "5 Yr" to 1., 5. etc.

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

    Output
    ------
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

    df = pd.read_html(str(table))[0]

    # Cols to keep
    maturity_cols = df.columns
    maturity_cols = maturity_cols[(maturity_cols.str.contains('Yr')) | maturity_cols.str.contains('Mo')]

    df = df[maturity_cols]

    renamer = {maturity_cols[i]: maturity_string_to_years(maturity_cols[i]) for i in range(len(maturity_cols))}
    df = df.rename(columns=renamer)

    # To generate a smoothed yield curve, take the average over each of the maturities and then regress
    Y = df.mean(axis=0) / 100

    curve, sts = calibrate_ns_ols(np.array(Y.index, dtype=float), Y.values, tau0=1.0)

    if sts.status > 0:
        Warning(f'Message: {sts.message}')

    # Curve object can be used to provide smoothness
    return curve
