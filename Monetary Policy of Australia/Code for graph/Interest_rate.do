/*

Ce document est utilisé pour produire un graphique traduisant les tendances
des années précédentes quant au taux d'intéret officiel de Reserve Bank of Australia
Nous voulons produire des données trimestrielle du taux d'intéret sur une base mensuelle

*/
clear

cd "C:\Users\Ribak\Documents\ProjetPersonnel\Monetary Policy of Australia\Data"

***retirer les lignes inutiles et vides dans la base de données***

import delimited "InterestRate", rowrange(266) colrange(1:3) clear





/* 
Affichage des données et clarification 
*/


label variable v1 "Date"
label variable v2 "Cash rate Target"
label variable v3 "Interbank Overnight Cash rate"


gen date = monthly(v1, "MY")

format date %tmCCYY

gen mois = month(date)

list v1 date v2 in 1/10

/* Graph monthly interest_rate */

twoway (line v2 date, color(orange)), ///
	title("Quarterly Cash Rate Target") ///
	ylabel(0(2)15, angle(0) labsize(small)) ///
	ytitle("Cash Rate Target (%)", size(medium)) ///
	xlabel(, angle(45)) ///
	xline(2008) ///
	note("1993 : Inflation Target introduced            Data source: Reserve Bank of    Australia", size(small) color(blue)) ///
	plotregion(margin(3 3 3 3)) ///  
	graphregion(margin(7 7 7 7)) /// 
	scheme(s1color)
 
graph export "Quaterly Interest rate monthly.pdf", replace

gen trimestre = quarter(date)
gen annee = year(date)

drop date
drop mois
drop trimestre
drop annee

/* Générer des trimestre pour avoir des données Trimestrielles, calculer la moyenne
des taux d'intéret pour chaque trimestre et en faire un graphique */


gen date = date(v1, "MY") 

format date %tm 

gen trimestre = quarter(date)
gen annee = year(date)


collapse (mean) interest_rate = v2, by(annee trimestre) ///Taux d'intéret par trimestre

format interest_rate %9.2f ///Forcer les 2 chiffres aprés la virgule

*list interest_rate

/*Formation du graphique*/



twoway (line interest_rate annee, sort lcolor(red)) /// 
       , title("Quarterly Cash Rate Target", size(large) lcolor(black)) ///  
         xtitle("Years", size(medium) lcolor(blue)) /// 
		 xline(1993, lcolor(green) lpattern(dash)) ///
         ytitle("Interest Rate (%)", size(medium) lcolor(blue)) ///  
         xlabel(1990(4)2024, angle(45) labsize(small)) ///  
         ylabel(0(2)15, angle(0) labsize(small)) ///  
         legend(off) ///  
         note("1993 : Inflation Target introduced            Data source: Reserve Bank of Australia", size(small) color(blue)) ///
         scheme(s1color) ///       
         xscale(range(1990 2024)) ///  
         bgcolor(white) ///  
         plotregion(margin(3 3 3 3)) ///  
         graphregion(margin(7 7 7 7)) ///  
         
graph export "Quaterly Interest rate 2.pdf", replace


