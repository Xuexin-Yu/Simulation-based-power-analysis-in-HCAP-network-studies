* -----------------------------------------------------------------------------
* Author: Xuexin Yu
* Date: September 30, 2025 (summarized on Oct 3 2026)
* Research project: Simulation-based power analyses in HCAP Network studies
* -----------------------------------------------------------------------------

* Table 1. Distribution of HCAP cognitive factor scores
use "HCAP_harmonized_data_PsyMCA-v6 (2).dta"
*Summarize cognitive factor score distribution 
sum fgcp  
by study: sum fgcp  

*Calculate age coefficients in each single HCAP study 
by study: reg fgcp age i.female 

*Calculate age coefficients in pooled analysis 
reg fgcp age i.female i.study 

* -----------------------------------------------------------------------------
*  Simulation code 
* -----------------------------------------------------------------------------

* Codebook
* npg: Sample size.
* sd: Standard deviation.
* delta: Effect size.
* proportion: Exposure prevalence

* -----------------------------------------------------------------------------
* Part I. Single-country analysis
* -----------------------------------------------------------------------------

* 1.1. HRS sample size = 3,347; effect size in HRS-HCAP SD units.
clear
set seed 3481
local reps = 10000

cap program drop sim1
program define sim1, rclass
    syntax, npg(string) delta(string) sd(string) proportion(string)
    drop _all
    local N = `npg'
    set obs `N'
    gen x = runiform() < `proportion'
    gen y = invnorm(uniform()) * `sd' + `delta' * x
    ttest y, by(x)
    return scalar p = `r(p)'
    return scalar hit = `r(p)' < 0.05
end

* Effect size = 0.05 SD

*Loop through exposure prevalence from 5%-50% in 5% increments 
foreach x of numlist 0.05(0.05)0.5 { 
    qui {
        simulate hit = r(hit), reps(`reps') nodots: ///
            sim1, npg(3347) delta(0.05) sd(1) proportion(`x')

        if `x' == 0.05 {
            noisily di _n _col(10) "HRS % exposed" _col(20) "Power"
        }

        su hit
        noisily di _col(10) `x' _col(20) %5.3f `r(mean)'
    }
}

* Change x, npg, delta, and sd for other scenarios.

* -----------------------------------------------------------------------------
* Part II. Additional sample size needed for power > 0.80
* -----------------------------------------------------------------------------

* Example study: HRS-HCAP.
clear
set seed 3481
local reps = 10000

cap program drop sim1
program define sim1, rclass
    syntax, npg(string) delta(string) sd(string) proportion(string)
    drop _all
    local N = `npg'
    set obs `N'
    gen x = runiform() < `proportion'
    gen y = invnorm(uniform()) * `sd' + `delta' * x
    ttest y, by(x)
    return scalar p = `r(p)'
    return scalar hit = `r(p)' < 0.05
end

* Effect size = 0.05 SD.

* Loop through sample size from 63347 to 65347 in 200 increments
foreach x of numlist 63347(200)65347 { 
    qui {
        simulate hit = r(hit), reps(`reps') nodots: ///
            sim1, npg(`x') delta(0.05) sd(1) proportion(0.05)

        if `x' == 3347 {
            noisily di _n _col(10) "HRS sample size" _col(20) "Power"
        }

        su hit
        noisily di _col(10) `x' _col(20) %5.3f `r(mean)'
    }
}

/* Tip: Restrict x to a reasonable range to reduce computation time. 
I originally started with a lower sample size with a larger increment to identify a reasonable range*/ 

* -------------------------------------------------------------------------------------
* Part III. Power calculation for effect heterogeneity in pooled cross-country analysis 
* -------------------------------------------------------------------------------------

* Example: Study 0 is CHARLS-HCAP, and study 1 is LASI-DAD.
clear
set seed 3481
local reps = 10000

* Effect heterogeneity = 0.05.
cap program drop sim
program define sim, rclass
    syntax, npg1(string) npg0(string) delta(string) sd1(string) sd0(string) ///
        delta2(string) delta1(string) prop1(string) prop0(string)
    drop _all
    local N = `npg1' + `npg0'
    set obs `N'

    * Generate the country indicator
    gen em = runiform() < (`npg1' / (`npg1' + `npg0'))

    * Set the treatment distribution within each level of the country indicator
    gen x = runiform() < `prop1' if em == 1
    replace x = runiform() < `prop0' if em == 0

    * Generate y based on the effect size within each country and the interaction term.
    gen y = invnorm(uniform()) * `sd0' + `delta' * x if em == 0
    replace y = invnorm(uniform()) * `sd1' + `delta' * x + ///
        `delta1' * em + `delta2' * x * em if em == 1

    * Set delta2 to represent effect heterogeneity.
    regress y x i.x##i.em
    matrix list r(table)
    local pval = r(table)[4, 9]
    return scalar pval = `pval'
    return scalar hit = `pval' < 0.05
end

*Loop through exposure prevalence in study 1 from 5-50% in 5% increments
*Exposure prevalence in study 0 fixed at 5% 
foreach x of numlist 0.05(0.05)0.5 {
    qui {
        simulate hit = r(hit), reps(`reps') nodots: ///
            sim, npg1(4096) npg0(9755) delta(0.0) delta1(-0.026) ///
            delta2(0.05) sd1(0.87) sd0(0.98) prop1(`x') prop0(0.05)

        su hit
        noisily di _col(10) `x' _col(20) %5.3f `r(mean)'
    }
}
* Change x, npg, delta, and sd for other scenarios.

* ------------------------------------------------------------------------------
* Part IV. Single-country analysis using internally standardized SD units
* ------------------------------------------------------------------------------

* Internal standardization means that cognitive factor z-scores have a mean of 0 and an SD of 1 in the model.

* 1.1. ELSA sample size = 1273
clear
set seed 3481
local reps = 10000

cap program drop sim1
program define sim1, rclass
    syntax, npg(string) delta(string) sd(string) proportion(string)
    drop _all
    local N = `npg'
    set obs `N'
    gen x = runiform() < `proportion'
    gen y = invnorm(uniform()) * `sd' + `delta' * x
    ttest y, by(x)
    return scalar p = `r(p)'
    return scalar hit = `r(p)' < 0.05
end

* Effect size = 0.05 SD.
foreach x of numlist 0.05(0.05)0.5 {
    qui {
        simulate hit = r(hit), reps(`reps') nodots: ///
            sim1, npg(1273) delta(0.05) sd(1.0) proportion(`x')

        if `x' == 0.05 {
            noisily di _n _col(10) "% exposed" _col(20) "Power"
        }

        su hit
        noisily di _col(10) `x' _col(20) %5.3f `r(mean)'
    }
}
* SD=1.07 for ELSA-HCAP in Part I. SD=1 here because cognitive factor scores were internally standardized. 

* -----------------------------------------------------------------------------
* Part V. Sample-size calculation using internally standardized SD units
* -----------------------------------------------------------------------------

* The same approach applies to all HCAP studies with a mean of 0 and an SD of 1.
clear
set seed 3481
local reps = 10000

cap program drop sim1
program define sim1, rclass
    syntax, npg(string) delta(string) sd(string) proportion(string)
    drop _all
    local N = `npg'
    set obs `N'
    gen x = runiform() < `proportion'
    gen y = invnorm(uniform()) * `sd' + `delta' * x
    ttest y, by(x)
    return scalar p = `r(p)'
    return scalar hit = `r(p)' < 0.05
end

* Effect size = 0.05 SD.
foreach x of numlist 100(100)1000 {
    qui {
        simulate hit = r(hit), reps(`reps') nodots: ///
            sim1, npg(`x') delta(0.05) sd(1) proportion(0.5)

        if `x' == 3347 {
            noisily di _n _col(10) "HRS sample size" _col(20) "Power"
        }

        su hit
        noisily di _col(10) `x' _col(20) %5.3f `r(mean)'
    }
}
*Tip: make x in a reasonable range to reduce compuation time

* -----------------------------------------------------------------------------
* Part VI. Cross-country analysis using internally standardized SD units
* -----------------------------------------------------------------------------

clear
set seed 3481
local reps = 10000

* Effect heterogeneity = 0.05.
cap program drop sim
program define sim, rclass
    syntax, npg1(string) npg0(string) delta(string) sd1(string) sd0(string) ///
        delta2(string) delta1(string) prop1(string) prop0(string)
    drop _all
    local N = `npg1' + `npg0'
    set obs `N'

    * Generate the country indicator
    gen em = runiform() < (`npg1' / (`npg1' + `npg0'))

    * Set the treatment distribution within each level of the country indicator
    gen x = runiform() < `prop1' if em == 1
    replace x = runiform() < `prop0' if em == 0

    * Generate y based on the effect size within each country and the interaction term.
    gen y = invnorm(uniform()) * `sd0' + `delta' * x if em == 0
    replace y = invnorm(uniform()) * `sd1' + `delta' * x + ///
        `delta1' * em + `delta2' * x * em if em == 1

    * Set delta2 to represent effect heterogeneity.
    regress y x i.x##i.em
    matrix list r(table)
    local pval = r(table)[4, 9]
    return scalar pval = `pval'
    return scalar hit = `pval' < 0.05
end

* Example: Effect difference = 0.05.
foreach x of numlist 0.05(0.05)0.5 {
    qui {
        simulate hit = r(hit), reps(`reps') nodots: ///
            sim, npg1(4096) npg0(9755) delta(0.0) delta1(-0.026) ///
            delta2(0.05) sd1(1) sd0(1) prop1(`x') prop0(0.05)

        su hit
        noisily di _col(10) `x' _col(20) %5.3f `r(mean)'
    }
}


