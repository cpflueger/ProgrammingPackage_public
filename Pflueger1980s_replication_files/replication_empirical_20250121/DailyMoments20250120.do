*"Back to the 1980s or Not?"
*cpflueger@uchicago.edu
*1/20/2025

*Generate Figures 9 and 10 and empirical moments in Table 3 

*set main path
global mainpath "C:\Users\cpflueger\Documents\GitHub\InflationExpectationsFormation\InflationExpectationsFormation\replication_empirical_20250121"

set scheme AC, permanently

****** Load GSW yields ******
* Import TIPS database 
import excel "$mainpath\feds200805_4.xlsx", sheet("feds200805") firstrow clear
drop BETA* TAU*

* Order database
gsort Date

* Destring variables
destring TIPSY02, replace   ignore("NA")
destring TIPSY03, replace   ignore("NA")
destring TIPSY04, replace   ignore("NA")
destring TIPSY05, replace   ignore("NA")
destring TIPSY06, replace   ignore("NA")
destring TIPSY07, replace   ignore("NA")
destring TIPSY08, replace   ignore("NA")
destring TIPSY09, replace   ignore("NA")
destring TIPSY10, replace   ignore("NA")
destring TIPSY19, replace   ignore("NA")
destring TIPSY20, replace   ignore("NA")

destring BKEVEN02, replace   ignore("NA")
destring BKEVEN03, replace   ignore("NA")
destring BKEVEN04, replace   ignore("NA")
destring BKEVEN05, replace   ignore("NA")
destring BKEVEN06, replace   ignore("NA")
destring BKEVEN07, replace   ignore("NA")
destring BKEVEN08, replace   ignore("NA")
destring BKEVEN09, replace   ignore("NA")
destring BKEVEN10, replace   ignore("NA")
destring BKEVEN19, replace   ignore("NA")
destring BKEVEN20, replace   ignore("NA")

* Save TIPS database
save "$mainpath\TIPS4.dta", replace

* Import Nominal Bond yields database 
import excel "$mainpath\feds200628_4.xlsx", sheet("Yields") firstrow clear
drop BETA* TAU*

* Order database
gsort Date

* Destring variables
destring SVENY01, replace   ignore("NA")
destring SVENY02, replace   ignore("NA")
destring SVENY03, replace   ignore("NA")
destring SVENY04, replace   ignore("NA")
destring SVENY05, replace   ignore("NA")
destring SVENY06, replace   ignore("NA")
destring SVENY07, replace   ignore("NA")
destring SVENY08, replace   ignore("NA")
destring SVENY09, replace   ignore("NA")
destring SVENY10, replace   ignore("NA")
destring SVENY19, replace   ignore("NA")
destring SVENY20, replace   ignore("NA")
destring SVENY30, replace   ignore("NA")

destring SVENF01, replace   ignore("NA")
destring SVENF02, replace   ignore("NA")
destring SVENF03, replace   ignore("NA")
destring SVENF04, replace   ignore("NA")
destring SVENF05, replace   ignore("NA")
destring SVENF06, replace   ignore("NA")
destring SVENF07, replace   ignore("NA")
destring SVENF08, replace   ignore("NA")
destring SVENF09, replace   ignore("NA")
destring SVENF10, replace   ignore("NA")
destring SVENF19, replace   ignore("NA")
destring SVENF20, replace   ignore("NA")
destring SVENF30, replace   ignore("NA")

* Save Nominal Bond yields database
save "$mainpath\Nominals4.dta", replace

****** Load S&P500 from St Louis Fred ******
import excel "$mainpath\SP500_3.xlsx", sheet("FRED Graph") firstrow clear
rename DATE Date
keep if SP500~=.
keep if SP500>1000
keep if SP500>0
* Order database 
gsort Date
tsset Date
* Generate S&P500 returns
gen retsp=(SP500-L.SP500)/L.SP500
* Save S&P500 database
save "$mainpath\SP500_3.dta", replace

****** Load Fed funds rate ******
import excel "$mainpath\FEDFUNDS.xls", sheet("FRED Graph") firstrow clear
lab var FEDFUNDS "Fed Funds"
gen ym=mofd(Date)
* Order database 
gsort ym 
tsset ym
* Save Fed funds rate database
save "$mainpath\FEDFUNDS.dta", replace

****** Load Inflation data from FRED (downloaded on 11/06/2024) ******
import excel "$mainpath\CPIAUCSL_updated.xls", sheet("FRED Graph") firstrow clear
lab var CPI "CPI Inflation (Ann. %)"
gen ym=mofd(Date)
* Order database 
gsort ym 
tsset ym
* Save Inflation database
save "$mainpath\CPI.dta", replace

******* Merge databases ***********
use "$mainpath\Nominals4.dta", replace
gsort Date
merge 1:1 Date using "$mainpath\TIPS4.dta"
drop _merge
merge 1:1 Date using "$mainpath\SP500_3.dta"
drop _merge
merge 1:1 Date using "$mainpath\CPI.dta"
drop _merge
merge 1:1 Date using "$mainpath\FEDFUNDS.dta"
drop _merge

* Complete Value-Weighted Return missing observatiosn with S&P500 returns
*replace vwretd=retsp if vwretd==.
rename retsp vwretd

* order database
rename Date date 
gsort date

* Gnerate year and quarter variables
gen year = year(date)
gen quarter = qofd(date)

* Save Final database
save "$mainpath/gsw_vwretd.dta", replace


************** Figures ******************
* Load Final database
use "$mainpath/gsw_vwretd.dta", clear

* Generate labels
lab var TIPSY10 "TIPS 10 YR"
lab var SVENY10 "Nominal 10 YR"
lab var BKEVEN10 "Breakeven 10 YR"
lab var BKEVEN05 "Breakeven 5 YR"
// lab var CPI "CPI Inflation"
lab var FEDFUNDS "Fed Funds Rate"
lab var date "Date"

*Figure 9

****************** Updated Figure postcovid *************************
graph twoway (line CPI FEDFUNDS date if year >= 2018 & year<2024, ///
    lwidth(thick thick) xline(21985, lcolor(black) lpattern(dash)) lpattern(solid dash) ///
    xsc(r(21185 23376)) ///  
	legend(pos(6) col(2) size(4)) ///
    xlabel(21185 21550 21915 22281 22646 23011 23376, format(%tdCCYY) nogrid) /// 
    xtitle("") ylab(,nogrid) lcolor(red blue))
graph export "$mainpath\postcovid_Pflueger.pdf", as(pdf) name("Graph") replace
*********************************************************************

*some additional plots of the raw data. These are not in the paper, but in some versions of the slides
graph twoway (line TIPSY10 date if year>=2018, yaxis(1) xline(21985) ) (line SP500 date if year>=2018, yaxis(2))
graph export "$mainpath\plot_levels_TIPS_postpandemic_Oct2023.pdf", as(pdf) name("Graph") replace

graph twoway (line SVENY10 date if year>=2018, yaxis(1) xline(21985) ) (line TIPSY10 date if year>=2018, yaxis(2))
graph export "$mainpath\plot_levels_nom_real_postpandemic_Oct2023.pdf", as(pdf) name("Graph") replace

graph twoway (line SVENY10 date if year>=2018, yaxis(1) xline(21985) ) (line SP500 date if year>=2018, yaxis(2))
graph export "$mainpath\plot_levels_nom_postpandemic_Oct2023.pdf", as(pdf) name("Graph") replace

graph twoway (line BKEVEN10 date if year>=2018, yaxis(1) xline(21985) ) (line SP500 date if year>=2018, yaxis(2))
graph export "$mainpath\plot_levels_postpandemic_Oct2023.pdf", as(pdf) name("Graph") replace

graph twoway (line BKEVEN10 BKEVEN05 date if year>=2018, yaxis(1) xline(21985) ) 
graph export "$mainpath\breakeven_only_Oct2023.pdf", as(pdf) name("Graph") replace
************ Generate term structure of covariances ************
* Rename variables for better consistency and readability
forvalues i=2/9 {
	rename SVENY0`i' SVENY`i'
	rename SVENF0`i' SVENF`i'
	rename TIPSY0`i' TIPSY`i'
	rename TIPSF0`i' TIPSF`i'
	rename BKEVEN0`i' BKEVEN`i'
	rename BKEVENF0`i' BKEVENF`i'
}

* Set data as tim-series on a daily frequency
tsset date

* Generate new variables for term structure calculations
forvalues i=2/20 {
	gen ry`i' = -s.SVENY`i'*`i'
	gen rtips`i' = -s.TIPSY`i'*`i'
	gen rf`i' = -s.SVENF`i'
	gen rby`i' = -s.BKEVEN`i'*`i'
	gen rbf`i' = -s.BKEVENF`i'
}

* Adjust market return variable to express it in percentage terms
replace vwretd=vwretd*100

forvalues i=2/20 {
	gen cov_ry`i'= ry`i'*vwretd
	gen cov_rf`i'= rf`i'*vwretd
	gen cov_by`i'= rby`i'*vwretd
	gen cov_bf`i'= rbf`i'*vwretd
	gen cov_bry`i'= rby`i'*ry`i'
	gen cov_brf`i'= rbf`i'*rf`i'
	gen var_ry`i' = ry`i'^2
	gen var_rf`i' = rf`i'^2
	gen var_by`i' = rby`i'^2
	gen var_bf`i' = rbf`i'^2
	gen cov_ryrf`i' = ry`i'*rf`i'
}

********* Compute daily bond yields term structure *********
* Load Final database
use "$mainpath/gsw_vwretd.dta", clear

* Rename variables for better consistency and readability
forvalues i=2/9 {
	rename SVENY0`i' SVENY`i'
	rename SVENF0`i' SVENF`i'
	rename TIPSY0`i' TIPSY`i'
	rename TIPSF0`i' TIPSF`i'
	rename BKEVEN0`i' BKEVEN`i'
	rename BKEVENF0`i' BKEVENF`i'
}

* Set data as tim-series on a daily frequency
tsset date

* Generate variables for term structure calculations
forvalues i=2/20 {
	gen ry`i' = -s.SVENY`i'*`i'
	gen rtips`i' = -s.TIPSY`i'*`i'
	gen rf`i' = -s.SVENF`i'
	gen rby`i' = -s.BKEVEN`i'*`i'
	gen rbf`i' = -s.BKEVENF`i'
}

// * Generate the inflation swap rate and scale it
// gen rswap=-10*s.InflationSwap10YR

* Scale market return variable for consistency (as a percentage)
replace vwretd=vwretd*100

* Generate month and quarter variables
gen month=mofd(date)
gen yq=qofd(date)

* Table 3
* Post-pandemic bond-stock betas using daily bond and stock returns 
* Early post-pandemic
regress ry10 vwretd if month>=723&month<=761
regress rtips10 vwretd if month>=723&month<=761

* Late post-pandemic
regress ry10 vwretd if month>761&year<=2023
regress rtips10 vwretd if month>761&year<=2023


* Keep only data from 2018 onward
keep if year>=2018

************** Rolling Regressions **************
* Loop through a list of variables for rolling regression and plot generation
foreach var of varlist ry10 rtips10 rby10 {

preserve
gsort date // Sort data by `date`
tsset date // Set `date` as the time-series variable

* Perform a rolling regression over a 120-day window
rolling _b _se, window(120) clear: regress `var' vwretd, ro

* Rename regression outputs for clarity
rename end date
rename _b_vwretd b_`var'
rename _se_vwretd se_`var'

* Create upper and lower bounds for a 90% confidence interval
gen b_`var'_upper=b_`var'+1.64*se_`var'
gen b_`var'_lower=b_`var'-1.64*se_`var'

* Label the variables for better graph interpretation
lab var b_`var' "Beta `var'"
lab var b_`var'_upper "90% CI"
lab var b_`var'_lower "90% CI"

* Sort and save the results
gsort date
save "$mainpath\rolling_`var'.dta", replace
restore
	
}

******* Merge databases ***********
* Merge results from different rolling regressions for combined analysis
gsort date 
merge date using "$mainpath\rolling_ry10.dta"
drop _merge
gsort date
merge date using "$mainpath\rolling_rtips10.dta"
drop _merge
gsort date
merge date using "$mainpath\rolling_rby10.dta"
drop _merge
gsort date

* Create lagged versions of variables to handle missing data in rolling analysis
gen b_rtips10_lag=L120.b_rtips10
gen b_ry10_lag=L120.b_ry10
gen b_rby10_lag=L120.b_rby10

* Replace values with missing if lagged values are missing 
replace b_rtips10=. if b_rtips10_lag==.
replace b_ry10=. if b_ry10_lag==.
replace b_rby10=. if b_rby10_lag==.

replace b_rtips10_upper=. if b_rtips10_lag==.
replace b_ry10_upper=. if b_ry10_lag==.
replace b_rby10_upper=. if b_rby10_lag==.

replace b_rtips10_lower=. if b_rtips10_lag==.
replace b_ry10_lower=. if b_ry10_lag==.
replace b_rby10_lower=. if b_rby10_lag==.

* Label final variables for interpretability in plots and analyses
lab var b_ry10 "Nominal Bond Beta"
lab var b_rtips10 "Infl-Indexed Bond Beta"
lab var b_rby10 "Breakeven-Stock Beta"
lab var date "Date"
gen zero=0 // Generate a zero line for plotting reference

lab var b_rtips10 "TIPS Bond-Stock Beta"
lab var b_ry10 "Nominal Bond-Stock Beta"
lab var b_rby10 "Breakeven-Stock Beta"
lab var b_ry10_upper "90% CI"
lab var b_ry10_lower "90% CI"

* Save the final combined dataset
save "$mainpath/betas_rolling.dta", replace


************** Figures ******************
* Load the final combined dataset
use "$mainpath/betas_rolling.dta", clear

gen date1=mdy(06,30,2022)
save "$mainpath\beta_nom_rolling_8.dta", replace
use "$mainpath\beta_nom_rolling_8.dta", clear
gen date2=mdy(06,30,2018)
keep if date>=date2

*Figure 10, Panel A
****************** Updated Figure plot_breakeven_postpandemic_Feb2024 ************************
graph twoway (line b_rby10 b_rby10_upper b_rby10_lower date if b_ry10 ~= . & year<2024, ///
    lcolor(blue black black green black black) ///
    lwidth(thick thin thin thick) ///
    lpattern(solid dash_dot dash_dot solid dash_dot dash_dot) ///
    xline(21985, lcolor(black)) ///
    yline(0, lcolor(black)) ///
    legend(pos(6) col(3) order(1 "Nominal-Minus-Real (Breakeven) Bond Beta" 2 "90% CI") size(4)) ///
    ylab(-0.3(0.1)0.2, nogrid) ///
    ysc(r(-0.3 0.2)) ///
    xscale(range(21185 23200)) ///  // Adjust these numbers to crop the x-axis and reduce white space
    xtitle("") xlabel(21185 21550 21915 22281 22646 23011 23376, format(%tdCCYY) nogrid) ///  // Ensure labels are within the cropped range
)
graph export "$mainpath\plot_breakeven_postpandemic_Feb2024_Pflueger.pdf", as(pdf) name("Graph") replace
***********************************************************************************************

*Figure 10, Panel B
****************** Updated Figure plot_betas_postpandemic_Feb2024 ************************
graph twoway (line b_ry10 b_rtips10 zero date if b_ry10 ~= . & year<2024, ///
    lcolor(red blue black black black black) ///
    lwidth(thick thick) ///
    lpattern(solid dash dash_dot dash_dot dash_dot dash_dot) ///
    xline(21985, lcolor(black)) ///
    xscale(range(21185 23200)) ///  // Adjust these numbers to crop the x-axis and reduce white space
    yline(0, lcolor(black)) ///
    legend(pos(6) col(3) order(1 "Nominal Bond-Stock Beta" 2 "Real Bond-Stock Beta") size(4)) ///
    yaxis(1) ///
    xlabel(21185 21550 21915 22281 22646 23011 23376, format(%tdCCYY) nogrid) ///  // Ensure labels are within the cropped range
    xtitle("") ylab(, nogrid) ///
)

graph export "$mainpath\plot_betas_postpandemic_Feb2024_Pflueger.pdf", as(pdf) name("Graph") replace
*******************************************************************************************
