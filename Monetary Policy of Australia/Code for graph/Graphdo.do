cd "C:\Users\Ribak\Documents\ProjetPersonnel\Monetary Policy of Australia\Data"


import excel "Data.xlsx", firstrow clear

format Year %tyCCYY

twoway(line Outputgap CashRate Inflation Year), ///
	title(" Output gap, Inflation and Interest Rate", margin(0 0 5 0)size(medium)) ///
	yline(0, lpattern(solid)) ///
	xline(2008, lcolor(gs11)) ///
	xlabel(1994(3)2024, angle(45) nogrid) ///
	xtitle("Year") ///
	ysize(300pt) ///
	ytitle("Rate in %") ///
	ylabel(-0.04 "-4%" -0.02 "-2%" 0 "0%" 0.02 "2%" 0.04 "4%" 0.06 "6%" 0.08 "8%" 0.10 "10%"  0.12 "12%" , ) ///
	bgcolor(red) ///
	graphregion(color(white)) ///
    plotregion(margin(0 0 0 0)color(white)) ///
	note("Source : IMF, ABS, RBA", size(vsmall) color(blue)) ///
	legend(position(below) order(1 "Output gap" 2 "Cash Rate" 3 "Inflation") row(1)) 
	
graph export "Inflation, Outputgap and Interest Rate.pdf", replace

/*************************************************
***création du graphique Cash rate and actual rate*** 
*************************************************/

*  Create graph prior 2008 (2008 included)

gen crisis = 2009
reg CashRate Inflation Outputgap CashRateLag if Year <= crisis
predict yhat

twoway line CashRate Year if Year <= crisis, lpattern(solid) ///
       || line yhat Year if Year <= crisis, lpattern(dash) ///
       title("Actual cash rate and Cash rate Target before 2008") ///
       xlabel(1994(3)2009, angle(45) nogrid) ///
       xtitle("Year") ///
       xsize(8) ///
       ytitle("Rate in %") ///
       ylabel(-0.01 "-1%" 0 "0%" 0.02 "2%" 0.04 "4%" 0.06 "6%" 0.08 "8%") ///
       yline(0, lpattern(solid)) ///
       legend(position(below) order(1 "Cash rate target" 2 "Actual Cash rate") row(1)) ///
       bgcolor(red) ///
       graphregion(color(white)) ///
       plotregion(color(white)) ///
       note("Source : IMF, ABS, RBA", size(small) color(blue))
	
graph export "Graph Actual rate and Rate Target before 2008.pdf", replace

* Create graph after 2008 (2008 excluded)

drop yhat
reg CashRate Inflation Outputgap CashRateLag if Year > crisis
predict yhat


twoway line CashRate Year if Year >= 2008, lpattern(solid) ///
       || line yhat Year if Year >= 2008, lpattern(dash) ///
       title("Actual cash rate and Cash rate Target after 2008") ///
       xlabel(2008(3)2024, angle(45) nogrid) ///
       xtitle("Year") ///
       xsize(8) ///
       ytitle("Rate in %") ///
       ylabel(-0.01 "-1%" 0 "0%" 0.02 "2%" 0.04 "4%" 0.06 "6%" 0.08 "8%") ///
       yline(0, lpattern(solid)) ///
       legend(position(below) order(1 "Cash rate target" 2 "Actual Cash rate") row(1)) ///
       bgcolor(red) ///
       graphregion(color(white)) ///
       plotregion(color(white)) ///
	   legend(position(below) order(1 "Cash rate target" 2 "Actual Cash rate") row(1)) ///
       note("Source : IMF, ABS, RBA", size(small) color(blue))
	   
graph export "Graph Actual rate and Rate Target After 2008.pdf", replace


/*    Estimating Taylor rule 1993       
 rt = r* + beta*(pi_t - pi_star) + gamma*xt 
 */


gen beta = 0.5
gen gamma = 0.5

gen inflation_target = 0.02
gen r_star = 0.02


gen r_taylor_1993 = r_star + (beta*(Inflation - inflation_target)) + gamma*Outputgap

*  Create overall graph including taylor rule 1993 ****

list r_taylor_1993

di "Estimation of taylor rule 1993 rate r* : " r_taylor_1993

twoway line CashRate Year , lpattern(solid) ///
       || line yhat Year , lpattern(dash) ///
	   || line r_taylor_1993 Year, lpattern(solid) ///
       title("Actual cash rate, Cash rate Target Taylor rule 1993") ///
       xlabel(1994(3)2024, angle(45) nogrid) ///
       xtitle("Year") ///
       xsize(8) ///
       ytitle("Rate in %") ///
       ylabel(-0.01 "-1%" 0 "0%" 0.02 "2%" 0.04 "4%" 0.06 "6%" 0.08 "8%") ///
       yline(0, lpattern(solid)) ///
       legend(position(below) order(1 "Cash rate target" 2 "Actual Cash rate") row(1)) ///
       bgcolor(red) ///
       graphregion(color(white)) ///
       plotregion(color(white)) ///
	   legend(position(below) order(1 "Cash rate target" 2 "Actual Cash rate" 3 " Taylor rule 1993 ") row(1)) ///
       note("Source : IMF, ABS, RBA", size(small) color(blue))




