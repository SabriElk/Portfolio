import yfinance as yf
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# Import Statsmodels
from statsmodels.tsa.api import VAR
from statsmodels.tsa.stattools import adfuller
from statsmodels.tools.eval_measures import rmse, aic
from statsmodels.tsa.stattools import grangercausalitytests

#This model is used to estimate the Value-at-Risk (VaR) of a portfolio of stocks using a Vector Autoregressive (VAR) approach.
#Download data from Yahoo Finance and compute returns
tickers = ['AAPL', 'MSFT', 'GOOGL', 'AMZN', 'TSLA','META']
rawdata = yf.download(tickers, start='2014-01-01', end='2024-01-01')
data = rawdata["Close"]
for i in tickers:
    Changes = data.pct_change().dropna()

#We use the return instead of price for the VAR model estimation because returns are typically stationary, while prices are not. Prices follows a random walk distribution.
#Stationarity is a key assumption for VAR models to ensure reliable and valid results.
#Test of causation using Granger causality test
#Granger causality tests whether one time series can predict another.

maxlag = 15
test = 'ssr_chi2test'
def grangers_causation_matrix(data, variables, maxlag, test='ssr_chi2test', verbose=False):    
    df = pd.DataFrame(np.zeros((len(variables), len(variables))), columns=variables, index=variables)
    for c in df.columns:
        for r in df.index:
            if r != c:
                test_result = grangercausalitytests(data[[r, c]], maxlag=maxlag, verbose=False)
                p_values = [round(test_result[i+1][0][test][1],4) for i in range(maxlag)]
                min_p_value = np.min(p_values)
                df.loc[r, c] = min_p_value
    df.columns = [var + '_x' for var in variables]
    df.index = [var + '_y' for var in variables]
    return df
print(grangers_causation_matrix(data=Changes, variables=Changes.columns, maxlag=maxlag, test=test))

#Since Amazon does not seem to cause any other stock, we will remove it from the dataset for the VAR estimation.
Changes_VAR = Changes.drop(columns=['AMZN'])

#Lets retest the causality matrix without Amazon
print(grangers_causation_matrix(data=Changes_VAR, variables=Changes_VAR.columns, maxlag=maxlag, test=test))