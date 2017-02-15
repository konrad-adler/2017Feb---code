******************************************************
* Monday, 28th April 2014
* Created by: Simon Fuchs
/* 

*/
******************************************************

************************************************************************************************************
*				 1. Prelims						   	 
************************************************************************************************************
clear all
*set memory 2g
*set matsize 800
set more off
global home "~/Dropbox/Research Shared/Movies/Analysis"
*global home C:/Users/Konrad/Dropbox/Movies/Analysis
cd "$home"

************************************************************************************************************
*				Create a balanced sample
************************************************************************************************************
use "$home/data_clean_Dec2016", replace
rename rev_dom us_rev
rename southaficaentireregion_rev safrica_rev

* create a balanced panel 1: select countries
keep fulltitle critics audience year production_budget sequel ///
australia_rev austria_rev belgium_rev brazil_rev  bulgaria_rev chile_rev   czechrepublic_rev  ///
finland_rev france_rev germany_rev greece_rev  hungary_rev iceland_rev italy_rev  lebanon_rev mexico_rev netherlands_rev ///
newzealand_rev norway_rev poland_rev portugal_rev singapore_rev  spain_rev  turkey_rev unitedarabemirates_rev /// 
unitedkingdom_rev us_rev ///
major genre_coded sequel remake adaptation // NOTE: these measures probably have to be improved

local countries  australia_rev austria_rev belgium_rev brazil_rev bulgaria_rev  chile_rev   czechrepublic_rev  ///
finland_rev france_rev germany_rev greece_rev  hungary_rev iceland_rev italy_rev  lebanon_rev mexico_rev netherlands_rev ///
newzealand_rev norway_rev poland_rev portugal_rev singapore_rev  spain_rev  turkey_rev unitedarabemirates_rev /// 
unitedkingdom_rev us_rev 


* additional filters
drop if production_budget ==.
drop if year<2000
drop if year>2014


* create a balanced panel 2: drop movies with missing 
foreach x in `countries'{
 drop if `x' ==.
}


local countries  australia_rev austria_rev belgium_rev brazil_rev bulgaria_rev  chile_rev   czechrepublic_rev  ///
finland_rev france_rev germany_rev greece_rev  hungary_rev iceland_rev italy_rev  lebanon_rev mexico_rev netherlands_rev ///
newzealand_rev norway_rev poland_rev portugal_rev singapore_rev  spain_rev  turkey_rev unitedarabemirates_rev /// 
unitedkingdom_rev us_rev 


* # ************************************************************
*Firstly generate country specific market size per genre/share
* # ************************************************************


foreach x in `countries'{
gen S_`x'=0
bys year: egen X_`x'=sum(`x') /* X_c: total revenues per year in country c*/
replace S_`x'= `x'/X_`x' /* S_mc the share of movie m in yearly revenue in country c*/
}


keep S_* year 
keep if year==2011
drop year


*hist year
*outsheet using "C:/Users/Konrad/Dropbox/Movies/Code - Struct 1216/data/balanced_panel.csv", replace
outsheet using "~/Dropbox/Research Shared/Movies/Analysis/balanced_panel.csv", replace
*clear all


************************************************************************************************************
*				Create a balanced sample
************************************************************************************************************
use "$home/data_clean_Dec2016", replace
rename rev_dom us_rev
rename southaficaentireregion_rev safrica_rev

* create a balanced panel 1: select countries
keep fulltitle critics audience year production_budget sequel ///
australia_rev austria_rev belgium_rev brazil_rev  bulgaria_rev chile_rev   czechrepublic_rev  ///
finland_rev france_rev germany_rev greece_rev  hungary_rev iceland_rev italy_rev  lebanon_rev mexico_rev netherlands_rev ///
newzealand_rev norway_rev poland_rev portugal_rev singapore_rev  spain_rev  turkey_rev unitedarabemirates_rev /// 
unitedkingdom_rev us_rev ///
major genre_coded sequel remake adaptation // NOTE: these measures probably have to be improved

local countries  australia_rev austria_rev belgium_rev brazil_rev bulgaria_rev  chile_rev   czechrepublic_rev  ///
finland_rev france_rev germany_rev greece_rev  hungary_rev iceland_rev italy_rev  lebanon_rev mexico_rev netherlands_rev ///
newzealand_rev norway_rev poland_rev portugal_rev singapore_rev  spain_rev  turkey_rev unitedarabemirates_rev /// 
unitedkingdom_rev us_rev 


* additional filters
drop if production_budget ==.
drop if year<2000
drop if year>2014


* create a balanced panel 2: drop movies with missing 
foreach x in `countries'{
 drop if `x' ==.
}


local countries  australia_rev austria_rev belgium_rev brazil_rev bulgaria_rev  chile_rev   czechrepublic_rev  ///
finland_rev france_rev germany_rev greece_rev  hungary_rev iceland_rev italy_rev  lebanon_rev mexico_rev netherlands_rev ///
newzealand_rev norway_rev poland_rev portugal_rev singapore_rev  spain_rev  turkey_rev unitedarabemirates_rev /// 
unitedkingdom_rev us_rev 


* # ************************************************************
*Firstly generate country specific market size per genre/share
* # ************************************************************


foreach x in `countries'{
gen S_`x'=0
bys year: egen X_`x'=sum(`x') /* X_c: total revenues per year in country c*/
replace S_`x'= `x'/X_`x' /* S_mc the share of movie m in yearly revenue in country c*/
}


keep S_* year production_budget
keep if year==2011
drop year
keep production_budget

*hist year
*outsheet using "C:/Users/Konrad/Dropbox/Movies/Code - Struct 1216/data/balanced_panel.csv", replace
outsheet using "~/Dropbox/Research Shared/Movies/Analysis/balanced_panel_budget.csv", replace
*clear all
