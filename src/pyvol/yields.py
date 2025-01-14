'''
Helper file to download latest data from US treasury (latest month)
'''

from datetime import datetime
import pandas as pd
import requests
from bs4 import BeautifulSoup

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

cols_to_keep = ['Date'] + list(maturity_cols)
df = df[cols_to_keep]


# Renaming and simplifying the df
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


renamer = {maturity_cols[i]: maturity_string_to_years(maturity_cols[i]) for i in range(len(maturity_cols))}
df = df.rename(columns=renamer)
print(df.head())

# To generate a smoothed yield curve, take the average over each of the maturities and then regress
Y = df.mean(axis=0)
