/*Code for banking assignment : Tayulor rule */

clear

cd "D:\Master_Economie\Master_1\Central_banking_assignement\All plot"

import excel "Data.xlsx", firstrow clear

format Year %tyCCYY


/* Insight on interest rate smoothing :Taylor rule 
Taylor rule */
/* Goal is to estimate the simple form of Taylor rule */


reg CashRate Inflation Outputgap CashRateLag   

drop in 31/31
/* Display value for beta_hat, gamma_hat and p_hat */

di _b[Inflation]
di _b[Outputgap]
di _b[CashRateLag]

/* Calculate value of beta, gamma and p */

gen beta_hat = _b[Inflation]
gen gamma_hat = _b[Outputgap]
gen p_hat = _b[CashRateLag]


gen beta = beta_hat / (1-_b[CashRateLag])
gen gamma  = gamma_hat / (1-_b[CashRateLag])
gen p = p_hat

/* display value of betta, gamma and p */


di "Adujusted beta : " beta
di "Adujusted gamma : " gamma
di "Adujusted p  : " p

/* 5. */
/* Calculate pi_star */


gen real_rate = CashRate - Inflation
summarize real_rate
gen rr_star =  r(mean)


gen pi_star = (((1-p_hat)*rr_star)-_b[_cons]) / ((1-p_hat)*(beta-1))
di "Pi_star : " pi_star

di 1-p
di beta - 1
di (1-p)*(beta -1)
di (1-p)*rr_star - _b[_cons]

/* calculate r_star */

gen r_star = _b[_cons] / (1-p_hat)
di "Natural interest rate " r_star

*twoway(line CashRate Inflation Year)

***************************************************************
/* Study if there is a difference prior or after 2008 crisis */
***************************************************************

clear

import excel "Data.xlsx", firstrow clear

format Year %tyCCYY

gen crisis = 2008

reg CashRate Inflation Outputgap CashRateLag if Year <= crisis /* reg pour 1994-2008 */

/* Display value for beta_hat, gamma_hat and p_hat */

di _b[Inflation]
di _b[Outputgap]
di _b[CashRateLag]

/* Calculate value of beta, gamma and p */

gen beta_hat = _b[Inflation]
gen gamma_hat = _b[Outputgap]
gen p_hat = _b[CashRateLag]


gen beta = beta_hat / (1-_b[CashRateLag])
gen gamma  = gamma_hat / (1-_b[CashRateLag])
gen p = p_hat

/* display value of betta, gamma and p */


di "Adujusted beta : " beta
di "Adujusted gamma : " gamma
di "Adujusted p  : " p

/* 5. */
/* Calculate pi_star */


gen real_rate = CashRate - Inflation
summarize real_rate
gen rr_star =  r(mean)


gen pi_star = (((1-p_hat)*rr_star)-_b[_cons]) / ((1-p_hat)*(beta-1))
di "Pi_star : " pi_star

di 1-p
di beta - 1
di (1-p)*(beta -1)
di (1-p)*rr_star - _b[_cons]

/* calculate r_star */

gen r_star = _b[_cons] / (1-p_hat)
di "Natural interest rate " r_star

**************************

clear

import excel "Data.xlsx", firstrow clear

format Year %tyCCYY
gen crisis = 2008

reg CashRate Inflation Outputgap CashRateLag if Year > crisis /* Reg pour 2009-2023 */

/* Display value for beta_hat, gamma_hat and p_hat */

di _b[Inflation]
di _b[Outputgap]
di _b[CashRateLag]

/* Calculate value of beta, gamma and p */

gen beta_hat = _b[Inflation]
gen gamma_hat = _b[Outputgap]
gen p_hat = _b[CashRateLag]


gen beta = beta_hat / (1-_b[CashRateLag])
gen gamma  = gamma_hat / (1-_b[CashRateLag])
gen p = p_hat

/* display value of betta, gamma and p */


di "Adujusted beta : " beta
di "Adujusted gamma : " gamma
di "Adujusted p  : " p

/* 5. */
/* Calculate pi_star */


gen real_rate = CashRate - Inflation
summarize real_rate
gen rr_star =  r(mean)


gen pi_star = (((1-p_hat)*rr_star)-_b[_cons]) / ((1-p_hat)*(beta-1))
di "Pi_star : " pi_star

di 1-p
di beta - 1
di (1-p)*(beta -1)
di (1-p)*rr_star - _b[_cons]

/* calculate r_star */

gen r_star = _b[_cons] / (1-p_hat)
di "Natural interest rate " r_star

***************************
/* création du tableau */
***************************

clear

import excel "Data.xlsx", firstrow clear

format Year %tyCCYY

*ssc install outreg2

gen crisis = 2008

reg CashRate Inflation Outputgap CashRateLag
est store eq1

reg CashRate Inflation Outputgap CashRateLag if Year > crisis /* Reg pour 2009-2023 */
est store eq2

reg CashRate Inflation Outputgap CashRateLag if Year <= crisis /* reg pour 1994-2008 */
est store eq3

set more off
outreg2 [eq*] using regress, e(all) replace see excel

