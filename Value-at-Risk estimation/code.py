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

#Téléchargement des données de Yahoo Finance
tickers = ['AAPL', 'MSFT', 'GOOGL', 'AMZN', 'TSLA','META']
rawdata = yf.download(tickers, start='2014-01-01', end='2024-01-01')
data = rawdata["Close"]
for i in tickers:
    Changes = data.pct_change().dropna()

#Analyse descriptive des actions
sns.lineplot(data=data,dashes=False)
plt.title('Stock Prices Over Time')
plt.xlabel('Date')
plt.ylabel('Price in USD')
plt.savefig(r'C:\Users\Ribak\Documents\ProjetPersonnel\Value-at-Risk estimation\Graphs\stock_prices.pdf',dpi=300)
plt.close()

#Analyse des rendements cumulatives des actions
cumulative_returns = (1 + Changes).cumprod() - 1
sns.lineplot(data=cumulative_returns, dashes=False)
plt.title('Cumulative Returns Over Time')
plt.xlabel('Date')
plt.ylabel('Cumulative Return in %')
plt.savefig(r'C:\Users\Ribak\Documents\ProjetPersonnel\Value-at-Risk estimation\Graphs\cumulative_returns.pdf',dpi=300)
plt.close()

#Attribuer un portfolio avec un poids égal pour chaque action
weight = np.array([1/len(tickers)]*len(tickers))
#On veut aussi un portfolie non vide
portfolio_value = 100000

#calcul du retour de notre portfolio 
returns = (Changes * weight).sum(axis=1)
print(returns)

#Calcul des retours sur différentes périodes 
period = [7,31,90,180,365]

for x in range(len(period)):
    days = period[x]
    range_returns = returns.rolling(window=days).sum()
    range_returns = range_returns.dropna()
    #print(range_returns)

    #Estimation du VaR
    confidence_level = 0.95 # 95% de confiance
    VaR = -np.percentile(Changes[i], (1 - confidence_level) * 100)

    #Graph de la distribuhtion des rendements
    range_returns_dollars = range_returns * portfolio_value
    sns.histplot(range_returns_dollars, bins=25, kde=True)
    plt.title(f'Distribution of {days}-day Portfolio Returns')
    plt.xlabel('Returns in Dollars')
    plt.ylabel('Frequency')
    plt.savefig(rf'C:\Users\Ribak\Documents\ProjetPersonnel\Value-at-Risk estimation\Graphs\portfolio_returns_{days}_days.pdf',dpi=500)
    plt.close()




