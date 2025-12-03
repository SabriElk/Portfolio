import yfinance as yf
import fredapi as fa
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

#This section is dedicated to data cleaning and preprocessing for Value-at-Risk estimation.
#Download data from Yahoo Finance
#Compute daily returns for each stock

tickers = ['AAPL', 'MSFT', 'GOOGL', 'AMZN', 'TSLA','META']
rawdata = yf.download(tickers, start='2014-01-01', end='2024-01-01')
data = rawdata["Close"]
for i in tickers:
    Changes = data.pct_change().dropna()

#download data from FRED
#Since data from FRED requires an API key, make sure to set it up properly.
Fred_api = pd.read_csv(r'C:\Users\Ribak\Documents\ProjetPersonnel\API.csv')
key = Fred_api.iloc[0, 0]
print(key)

fred=fa.Fred(api_key=key)
#Download the 10-Year Treasury Constant Maturity Rate (DGS10) from FRED
treasury_data = fred.get_series('DGS10', observation_start='2014-01-01', observation_end='2024-01-01')
treasury_data = treasury_data.dropna()

#Plot the Treasury yield data
sns.lineplot(data=treasury_data)
plt.title('10-Year Treasury Constant Maturity Rate Over Time')
plt.xlabel('Date')
plt.ylabel('Yield in %')
plt.savefig(r'C:\Users\Ribak\Documents\ProjetPersonnel\Value-at-Risk estimation\Graphs\treasury_yield.pdf',dpi=300)
plt.close()