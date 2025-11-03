

cd "C:\Users\Ribak\Documents\ProjetPersonnel\Monetary Policy of Australia\Data"

import excel "outputgap", sheet(Data) clear

drop C

list A


label variable A "Date"
label variable B "Output gap in %"

/* création du graph */

twoway(line B A, color(red)), ///
	title("Output gap in %", margin(medium)) ///
	xlabel(1980(10)2024, angle(45) nogrid) ///
	xline(2008, lcolor(gs14) lpattern(solid)) ///
	text(2 2008 "    2008 Crisis     ", size(small)) ///
	xline(2018, lcolor(gs14) lpattern(solid)) ///
	text(2 2018 "COVID-19", lpattern(solid) size(small)) ///
	bgcolor(none) ///
	ylabel(, nogrid) ///
	yline(0, lcolor(black9) lpattern(solid)) ///
	note("Source : FMI", color(blue)) 
	
help xline