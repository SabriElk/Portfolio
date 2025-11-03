/* change directory */

cd "C:\Users\Ribak\Documents\ProjetPersonnel\Monetary Policy of Australia\Data"

import excel "Inflation.xlsx", sheet(Data1) cellrange(A12:AB315) clear

format A %td
/* AB est le variable qui correspond au changement CPI national */

label variable A "Date"
label variable AB "Percentage Change - CPI"
label data "Percentage change from Previous Period - CPI"

/* creation du graphique */

twoway (line AB A if A >= td(01jan1990) , color(red)), ///
	title("Percentage Change from Previous Period National - CPI ") ///
	xtitle(" Date") ///
	xlabel(,angle(45) format(%tdCCYY)) ///
	ytitle("Inflation rate (%)") ///
	yline(0) ///
	note("Source : Australian Bureau of Statistics", size(small) color(blue))
	
graph export "Percentage change - CPI.pdf"
