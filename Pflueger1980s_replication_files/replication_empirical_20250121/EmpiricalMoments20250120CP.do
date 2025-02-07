* empirical calibration moments for "Back to the 1980s or Not?"
*cpflueger@uchicago.edu
***
graph drop _all
global mainpath "C:\Users\cpflueger\Documents\GitHub\InflationExpectationsFormation\InflationExpectationsFormation\replication_empirical_20250121"
set scheme AC, permanently
import excel "$mainpath\DataTable_2022.xlsx", sheet("Data") firstrow clear

gen quarter=qofd(qdate)
gen year=yofd(qdate)
gsort quarter
format quarter %tq
tsset quarter

* Compute annualized log inflation
*GDP deflator used in local projections
gen inflq=400*(log(gdpdef)-log(L.gdpdef))
*CPI used to compute approximate TIPS returns
gen inflq_cpi=400*(log(CPI)-log(L.CPI))
lab var inflq "Quarterly Log Inflation (Ann %), Measured with GDP Deflator"

gen inflwage=400*(log(wage)-log(L.wage))
lab var inflwage "Quarterly Log Wage Inflation (Ann %)"
gen inflann=0.25*(inflq+L.inflq+L2.inflq+L3.inflq)
*CG forecast revision
gen forecast_revision = dpgdp4-L.dpgdp5
gen forecast_error = inflq-L3.dpgdp4

*5-year inflation
tssmooth ma infl5yr=inflq, window(20 0 0)
gen Delta_inflann=inflann-L4.inflann

*Log potential output
gen potx=100*log(gdppot)
*Log total output
gen totalx=100*log(gdp)

*log real consumption
gen c=log(Consumption)

*Log output gap relative to BAE potential output. This is in natural units and does not need to be scaled. 
*Percentage deviation from natural output
gen x=100*(log(gdp)-log(gdppot))

*******************************************************************************
*Sub-period cut offs
*******************************************************************************
gen q=qofd(qdate)

gen q1=78
gen q2= 165
gen q3= 208
*gen qgreenspan=110

drop if qdate==.
gsort qdate

gsort q
tsset q
***********************************************
*stochastically de-trended real output 
************************************************

tssmooth exponential y_exp3=log(gdp), parms(0.01) samp0(20)
gen y_se3=100*(log(gdp)-L.y_exp3)

tssmooth exponential y_exp2=log(gdp), parms(0.07) samp0(20)
gen y_se2=100*(log(gdp)-L.y_exp2)


lab var y_se2 "Detrended output phi=0.93"
lab var y_se3 "Detrended output phi=0.99"

*correlation between output gap and different versions of stochastically de-trended real output (phi=0.99 vs. phi=0.93)
pwcorr x y_se2 y_se3 if q>=q1&year<=2019

*volatility of Gilchrist-Zakrajsek spread in both subperiods (footnote 14)
sum gz_spread if q>=q1&q<q2
sum gz_spread if q>=q2&year<=2019

***************************************************
*Asset Pricing variables
***************************************************
*quarterly log nominal bond excess return
gen xr_bond=(-9.75*SVENY10+10*L.SVENY10-L.tbill/4)
*1-year bond excess return for Campbell-Shiller regressions
gen xr_bond1yr=(-9*SVENY09+10*L4.SVENY10-(L.tbill-L2.tbill-L3.tbill-L4.tbill)/4)
*Slope of the term structure
gen slope=SVENY10-tbill

*quarterly log equity excess return
gen xr_equity=(vwret-L.tbill)/4
*annual log equity excess return
gen xr_equity1yr=(xr_equity+L.xr_equity+L2.xr_equity+L3.xr_equity)/100

*US and UK log returns on indexed bonds
gen xr_TIPS=(-39*TIPSY10+40*L.TIPSY10-L.tbill+inflq_cpi)/4
gen xr_TIPSUK=(-39*linkeruk10+40*L.linkeruk10-L.tbill+inflq_cpi)/4

*log price-dividend ratio
replace pd=log(pd)
gen dp=-pd
lab var dp "Log Dividend-Price Ratio"
lab var pd "Log Price-Dividend Ratio"

lab var xr_equity "Quarterly Equity Excess Return"
lab var xr_bond "Quarterly 10YR Bond Excess Return"

***********************************************
*macro variables - annual change in Fed Funds rate
***********************************************
gen Deltai=FedFunds-L4.FedFunds
gen Deltac=100*(c-L4.c)
gen Deltainflq=inflq-L4.inflq
gen DeltaE10_CPI=SPF_CPI10-L4.SPF_CPI10

**************************************************
*summary stats for asset pricing moments
**************************************************

******************************************************
*Stocks
*********************************************************
*average equity premium
*multiply to obtain equity premium in annualized percent and add 0.5 variance for Jensen's inequality adjustment
gen xr_equity_scaled=(xr_equity+0.5*xr_equity^2/100)*4
egen equity_premium1=mean(xr_equity_scaled) if q>=q1&q<q2
egen equity_premium2=mean(xr_equity_scaled) if q>=q2&year<=2019
drop xr_equity_scaled

*equity volatility
*multiply by 2 to display std in annualized percent
gen xr_equity_scaled=xr_equity*2
egen equity_vol1=sd(xr_equity_scaled) if q>=q1&q<q2
egen equity_vol2=sd(xr_equity_scaled) if q>=q2&year<=2019
drop xr_equity_scaled

*Sharpe ratio
gen sharpe1=equity_premium1/equity_vol1 if q>=q1&q<q2
gen sharpe2=equity_premium2/equity_vol2 if q>=q2&year<=2019

sum equity_premium1 equity_premium2
sum equity_vol1 equity_vol2
sum sharpe1 sharpe2
drop equity_premium* equity_vol* sharpe*

*AR(1) coefficient dividend yield
regress pd L.pd if q>=q1&q<q2
regress pd L.pd if q>=q2&year<=2019

*1-year log excess stock return onto lagged pd
regress xr_equity1yr L4.pd if q>=q1&q<q2, robust
regress xr_equity1yr L4.pd if q>=q2&year<=2019, robust

****************************************************************
*Bonds
****************************************************************

*bond return volatility
*multiply by 2 to display std in annualized percent
gen xr_bond_scaled=xr_bond*2
egen bond_vol1=sd(xr_bond_scaled) if q>=q1&q<q2
egen bond_vol2=sd(xr_bond_scaled) if q>=q2&year<=2019
drop xr_bond_scaled
sum bond_vol1 bond_vol2

drop bond_vol* 

// sum Delta_SVENY10 if q>=q1&q<q2
// sum Delta_SVENY10 if q>=q2&year<=2019

*Campbell-Shiller bond return predictability regressions
*1-yr excess returns on log yield slope
*1-year excess returns in percent
newey xr_bond1yr L4.slope if q>=q1&q<q2, lag(4)
newey xr_bond1yr L4.slope if q>=q2&year<=2019, lag(4)

**********************************
*bond return betas
**********************************
regress xr_bond xr_equity if q>=q1&q<q2, ro
regress xr_bond xr_equity if q>=q2&year<=2019, ro

regress xr_TIPSUK xr_equity if q>=q1&q<q2, ro
regress xr_TIPS xr_equity if q>=q2&year<=2019, ro
**********************************************
*macro calibration statistics 
************************************************
sum Deltac if q>=q1&q<q2
sum Deltai if q>=q1&q<q2

sum Deltac if q>=q2&year<=2019
sum Deltai if q>=q2&year<=2019

*volatility long-term inflation expectations 
gen SPF_CPI10_combined=SPF_CPI10
replace SPF_CPI10_combined=BC_CPI10 if SPF_CPI10==.

gen DeltaE10_CPI_combined=SPF_CPI10_combined-L4.SPF_CPI10_combined
sum DeltaE10_CPI_combined if q>=q1&q<q2
sum DeltaE10_CPI_combined if q>=q2&year<=2019
**************************************************
*test for rationality of inflation expectations: Appendix Table A3
******************************************************
newey forecast_error L3.forecast_revision if q<q2, lag(4)
estimates store rationality_period2
newey forecast_error L3.forecast_revision if q>=q1&q<q2, lag(4)
estimates store rationality_period4
newey forecast_error L3.forecast_revision if q>=q2&year<=2019, lag(4)
estimates store rationality_period5
regress forecast_error L3.forecast_revision if q<q2, ro
estimates store rationality_period2_rsq
regress forecast_error L3.forecast_revision if q>=q1&q<q2, ro
estimates store rationality_period4_rsq
regress forecast_error L3.forecast_revision if q>=q2&year<=2019, ro
estimates store rationality_period5_rsq

esttab rationality_period2 rationality_period4 rationality_period5 rationality_period2_rsq rationality_period4_rsq rationality_period5_rsq using "$mainpath\output\rationality_regressions.csv", star(* .10 ** .05 *** .01)  se(a2) nomtitle replace r2(2)

*************************************************************
*empirical moments for Figure 2
*************************************************************
*gen Dinflq=D.inflq

* Choose impulse response horizon
local hmax = 20

/* Generate LHS variables for the LPs */
* levels
gen Dc=D.c
*end-of-quarter timing convention
gen Fx=F.x

forvalues h = 0/`hmax' {
	gen infl_`h' = f`h'.inflq
	*end-of-quarter timing convention
	gen x_`h' = f`h'.Fx
	gen FedFunds_`h' = f`h'.FedFunds
	gen inflq_`h' = f`h'.inflq
}

graph drop _all
cap drop all

* macro.moments.b' in Calibrationgridsearch_supply2_init.mat
* Choose impulse response horizon
local hmax = 20

eststo clear
cap drop b u d Quarters Zero
gen Quarters = _n-1 if _n<=`hmax'
gen Zero =  0    if _n<=`hmax'
gen b=0
gen u=0
gen d=0
forv h = 0/`hmax' {
	* levels
	 newey x_`h' inflq L.inflq if q>=q1&q<q2, lag(`h')
replace b = _b[inflq]                    if _n == `h'+1
replace u = _b[inflq] + 1.96* _se[inflq]  if _n == `h'+1
replace d = _b[inflq] - 1.96* _se[inflq]  if _n == `h'+1
eststo
}

twoway ///
(rarea u d  Quarters,  ///
fcolor(gs13) lcolor(gs13) lw(none) lpattern(solid)) ///
(line b Quarters, lcolor(blue) ///
lpattern(solid) lwidth(thick)) ///
(line Zero Quarters, lcolor(black)), legend(off) ///
ytitle("Percent", size(medsmall)) xtitle("Quarter", size(medsmall)) ///
graphregion(color(white)) plotregion(color(white))
graph export "$mainpath/IRF_x_pi_pre2001.pdf", as(pdf)  replace

* macro.moments.b' in Calibrationgridsearch_demand2_init.mat
* Choose impulse response horizon
local hmax = 20

eststo clear
cap drop b u d Quarters Zero
gen Quarters = _n-1 if _n<=`hmax'
gen Zero =  0    if _n<=`hmax'
gen b=0
gen u=0
gen d=0
forv h = 0/`hmax' {
	* levels
	 newey x_`h' inflq L.inflq if q>=q2&year<=2019, lag(`h')
replace b = _b[inflq]                    if _n == `h'+1
replace u = _b[inflq] + 1.96* _se[inflq]  if _n == `h'+1
replace d = _b[inflq] - 1.96* _se[inflq]  if _n == `h'+1
eststo
}

twoway ///
(rarea u d  Quarters,  ///
fcolor(gs13) lcolor(gs13) lw(none) lpattern(solid)) ///
(line b Quarters, lcolor(blue) ///
lpattern(solid) lwidth(thick)) ///
(line Zero Quarters, lcolor(black)), legend(off) ///
ytitle("Percent", size(medsmall)) xtitle("Quarter", size(medsmall)) ///
graphregion(color(white)) plotregion(color(white))
graph export "$mainpath/IRF_x_pi_post2001.pdf", as(pdf)  replace

*macro.moments.bi in Calibrationgridsearch_supply2_init.mat
* Choose impulse response horizon
local hmax = 20

eststo clear
cap drop b u d Quarters Zero
gen Quarters = _n-1 if _n<=`hmax'
gen Zero =  0    if _n<=`hmax'
gen b=0
gen u=0
gen d=0
forv h = 0/`hmax' {
	* levels
	 newey FedFunds_`h' inflq L.inflq if q>=q1&q<q2, lag(`h')
replace b = _b[inflq]                    if _n == `h'+1
replace u = _b[inflq] + 1.96* _se[inflq]  if _n == `h'+1
replace d = _b[inflq] - 1.96* _se[inflq]  if _n == `h'+1
eststo
}

twoway ///
(rarea u d  Quarters,  ///
fcolor(gs13) lcolor(gs13) lw(none) lpattern(solid)) ///
(line b Quarters, lcolor(blue) ///
lpattern(solid) lwidth(thick)) ///
(line Zero Quarters, lcolor(black)), legend(off) ///
ytitle("Percent", size(medsmall)) xtitle("Quarter", size(medsmall)) ///
graphregion(color(white)) plotregion(color(white))
graph export "$mainpath/IRF_i_pi_pre2001.pdf", as(pdf)  replace

*macro.moments.bi in Calibrationgridsearch_demand2_init.mat
* Choose impulse response horizon
local hmax = 20

eststo clear
cap drop b u d Quarters Zero
gen Quarters = _n-1 if _n<=`hmax'
gen Zero =  0    if _n<=`hmax'
gen b=0
gen u=0
gen d=0
forv h = 0/`hmax' {
	* levels
	 newey FedFunds_`h' inflq L.inflq if q>=q2&year<=2019, lag(`h')
replace b = _b[inflq]                    if _n == `h'+1
replace u = _b[inflq] + 1.96* _se[inflq]  if _n == `h'+1
replace d = _b[inflq] - 1.96* _se[inflq]  if _n == `h'+1
eststo
}

twoway ///
(rarea u d  Quarters,  ///
fcolor(gs13) lcolor(gs13) lw(none) lpattern(solid)) ///
(line b Quarters, lcolor(blue) ///
lpattern(solid) lwidth(thick)) ///
(line Zero Quarters, lcolor(black)), legend(off) ///
ytitle("Percent", size(medsmall)) xtitle("Quarter", size(medsmall)) ///
graphregion(color(white)) plotregion(color(white))
graph export "$mainpath/IRF_i_pi_post2001.pdf", as(pdf)  replace

*macro.moments.bix in Calibrationgridsearch_supply2_init.mat
* Choose impulse response horizon
local hmax = 20

eststo clear
cap drop b u d Quarters Zero
gen Quarters = _n-1 if _n<=`hmax'
gen Zero =  0    if _n<=`hmax'
gen b=0
gen u=0
gen d=0
forv h = 0/`hmax' {
	* levels
	 newey x_`h' FedFunds L.FedFunds if q>=q1&q<q2, lag(`h')
replace b = _b[FedFunds]                    if _n == `h'+1
replace u = _b[FedFunds] + 1.96* _se[FedFunds]  if _n == `h'+1
replace d = _b[FedFunds] - 1.96* _se[FedFunds]  if _n == `h'+1
eststo
}

twoway ///
(rarea u d  Quarters,  ///
fcolor(gs13) lcolor(gs13) lw(none) lpattern(solid)) ///
(line b Quarters, lcolor(blue) ///
lpattern(solid) lwidth(thick)) ///
(line Zero Quarters, lcolor(black)), legend(off) ///
ytitle("Percent", size(medsmall)) xtitle("Quarter", size(medsmall)) ///
graphregion(color(white)) plotregion(color(white))
graph export "$mainpath/IRF_x_i_pre2001.pdf", as(pdf)  replace


*macro.moments.bix in Calibrationgridsearch_demand2_init.mat
* Choose impulse response horizon
local hmax = 20

eststo clear
cap drop b u d Quarters Zero
gen Quarters = _n-1 if _n<=`hmax'
gen Zero =  0    if _n<=`hmax'
gen b=0
gen u=0
gen d=0
forv h = 0/`hmax' {
	* levels
	 newey x_`h' FedFunds L.FedFunds if q>=q2&year<=2019, lag(`h')
replace b = _b[FedFunds]                    if _n == `h'+1
replace u = _b[FedFunds] + 1.96* _se[FedFunds]  if _n == `h'+1
replace d = _b[FedFunds] - 1.96* _se[FedFunds]  if _n == `h'+1
eststo
}

twoway ///
(rarea u d  Quarters,  ///
fcolor(gs13) lcolor(gs13) lw(none) lpattern(solid)) ///
(line b Quarters, lcolor(blue) ///
lpattern(solid) lwidth(thick)) ///
(line Zero Quarters, lcolor(black)), legend(off) ///
ytitle("Percent", size(medsmall)) xtitle("Quarter", size(medsmall)) ///
graphregion(color(white)) plotregion(color(white))
graph export "$mainpath/IRF_x_i_post2001.pdf", as(pdf)  replace

*macro.moments.bwage in Calibrationgridsearch_demand2_init.mat
*output gap and wage inflation
local hmax = 20

eststo clear
cap drop b u d Quarters Zero
gen Quarters = _n-1 if _n<=`hmax'
gen Zero =  0    if _n<=`hmax'
gen b=0
gen u=0
gen d=0
forv h = 0/`hmax' {
	* levels
	 newey x_`h' inflwage L.inflwage if q>=q2&year<=2019, lag(`h')
replace b = _b[inflwage]                    if _n == `h'+1
replace u = _b[inflwage] + 1.96* _se[inflwage]  if _n == `h'+1
replace d = _b[inflwage] - 1.96* _se[inflwage]  if _n == `h'+1
eststo
}

twoway ///
(rarea u d  Quarters,  ///
fcolor(gs13) lcolor(gs13) lw(none) lpattern(solid)) ///
(line b Quarters, lcolor(blue) ///
lpattern(solid) lwidth(thick)) ///
(line Zero Quarters, lcolor(black)), legend(off) ///
ytitle("Percent", size(medsmall)) xtitle("Quarter", size(medsmall)) ///
graphregion(color(white)) plotregion(color(white))
graph export "$mainpath/IRF_x_wage_post2001.pdf", as(pdf)  replace

save "$mainpath/macro_data.dta", replace

* Figure 1
use "$mainpath/macro_data.dta", clear
preserve
gsort quarter
tsset quarter

* Rolling regressions for nominal bond betas
rolling _b _se, window(24): regress xr_bond xr_equity, ro
rename end quarter
rename _b_xr_equity b_bond
rename _se_xr_equity se_bond
lab var b_bond "Nominal Bond Beta"
keep quarter b_bond se_bond
gsort quarter
save "$mainpath/rolling_nominalbeta.dta", replace
restore

preserve
gsort quarter
tsset quarter

* Rolling regression for real bond betas (TIPS)
rolling _b _se, window(24): regress xr_TIPS xr_equity, ro
rename end quarter
rename _b_xr_equity b_TIPS
rename _se_xr_equity se_TIPS
lab var b_TIPS "Infl-Indexed Bond Beta"
keep quarter b_TIPS se_TIPS
gsort quarter
save "$mainpath/rolling_TIPSbeta.dta", replace
restore


preserve
gsort quarter
tsset quarter

* Rolling regression for UK real bond betas (UK TIPS)
rolling _b _se, window(24): regress xr_TIPSUK xr_equity, ro
rename end quarter
rename _b_xr_equity b_TIPSUK
rename _se_xr_equity se_TIPSUK
lab var b_TIPSUK "UK Linker Bond Beta"
keep quarter b_TIPSUK se_TIPSUK
gsort quarter
save "$mainpath/rolling_TIPSUKbeta.dta", replace
restore

* Merge rolling regression results into one dataset
use "$mainpath/macro_data.dta", clear
gsort quarter
merge quarter using  "$mainpath/rolling_nominalbeta.dta"
drop _merge
gsort quarter
merge quarter using  "$mainpath/rolling_TIPSbeta.dta"
drop _merge
gsort quarter
merge quarter using  "$mainpath/rolling_TIPSUKbeta.dta"
drop _merge
gsort quarter
tsset quarter

* Require at least 10 observations
gen L10b_bond=L10.b_bond
gen L10b_TIPS=L10.b_TIPS
gen L10b_TIPSUK=L10.b_TIPSUK
replace b_bond=. if L10b_bond==.
replace b_TIPS=. if L10b_TIPS==.
replace b_TIPSUK=. if L10b_TIPSUK==.
replace b_TIPS=b_TIPSUK if b_TIPS==.

* Format quarter variable and label beta variables
format q %tq
lab var q "Quarter"
lab var b_TIPS "Infl-Indexed Bond Beta"

gen zeros = 0 // Helper variable for visualization

**************************** Figure 1 *************************************************************
graph twoway (line b_bond b_TIPS q if q >= q1, ///
    lcolor(red blue black black black black) ///
    lwidth(thick thick) ///
    lpattern(solid dash dash_dot dash_dot dash_dot dash_dot) ///
    xline(240, lcolor(black) lpattern(dash)) xline(165, lcolor(black) lpattern(dash)) ///
    yline(0, lcolor(black) lpattern(dash)) ///
    legend(pos(6) cols(2) size(4)) ///
    xtitle("") ylab(, nogrid) xlabel(, format(%tqCCYY) nogrid) ///  // Format to display only the year for quarterly data
)
graph export "$mainpath/BondBetas_rolling_Dec2023_noCI_Pflueger.pdf", as(pdf) replace
***************************************************************************************************
