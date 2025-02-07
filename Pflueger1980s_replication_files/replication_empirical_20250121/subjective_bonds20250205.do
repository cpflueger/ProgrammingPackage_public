*Back to the 1980s or Not?
*cpflueger@uchicago.edu
*02/05/2025

*This file generates expected excess bond returns (data) for Table 2
global mainpath "C:\Users\cpflueger\Dropbox\Documents\Econ\Ideas\InflationSources\data\replication_empirical_20250121"

*Load consensus 12-month forecast of the 10-year Treasury rate from BC data
*Blue Chip Financial Forecast data is proprietary
*pseudo data set
import delimited "$mainpath\us_analyst_estimates_pseudo.csv", clear
*uncomment the next line for actual data set
*import delimited "$mainpath\us_analyst_estimates.csv", clear
gen date_var = monthly(ym, "YM")
format date_var %tm
drop ym
rename date_var month
gsort month
keep month con12
order month
* Save dataset
save "$mainpath\bc_timeseries.dta", replace



*Load GSW yields - these have par yields, which is convenient
import delimited "$mainpath\feds200628.csv", clear 
gen dd= date(date,"MDY")
drop date
rename dd date
keep if svenpy05~=.
gen month=mofd(date)
gsort month
format month %tm
by month: keep if _n==1
gsort month
isid month
save "$mainpath\GSW_monthly.dta", replace


use "$mainpath\GSW_monthly.dta", clear

merge 1:1 month using "$mainpath\bc_timeseries.dta"
sum _merge
drop _merge
gsort month

gen year=yofd(dofm(month))
gsort month
tsset month

*one-year subjective risk premium
gen dur11=(1-(1+svenpy11/100)^(-11))/(1-(1+svenpy11/100)^(-1))
gen xr_nom10yr_survey=dur11*svenpy11-(dur11-1)*f.con12-sveny01

*quarterly summary statistics 
gen yq=qofd(dofm(month))
gsort yq -month
by yq: keep if _n==1
*subjective bond excess return for period 1 - Table 2 in the paper
sum xr_nom10yr_survey if yq<165
*subjective bond excess return for period 2 - Table 2 in the paper
sum xr_nom10yr_survey if yq>=165&year<=2019