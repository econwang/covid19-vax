**M ===============================================================================================
**M # Stata code for the vaccination analysis with the CHARLS data
**M     prepared for the paper on Nature Medicine (2023/1/24)
**M     1. Vax data initiation
**M     2. Data appending with weights and health-related variables
**M     3. Description
**M     4. Regression
**M ===============================================================================================

// based on the weight variable INDV_weight_ad2, further adjusted for COVID non-response and post-stratifications
gl weightvar COVID_weight_post2 

gl fdata_21 // path of CHARLS module data
gl fdata_22 // path of CHARLS module data
gl cdata_22 // path of CHARLS module data
gl processed_data   // path of processed data
gl varset_data // path of processed data
gl working_data // path of saving temporary data

** =================================================================
** 1. Vax data initiation for the 2022 survey
** =================================================================

use $fdata_22/COVID, clear

// use remark data to adjust the dose required for completing the primary series
// the data related to individual privacy are not released, including remark data
sort n006_remark
gen maxdoses3 = 1 if ~mi(n006_remark)
replace maxdoses3 = 0 if n006_remark=="第三针是加强针"
replace maxdoses3 = 0 if n006_remark=="第三针为加强针"
replace maxdoses3 = 3 if maxdoses3==1 & n006_remark=="3针为基础针"
replace maxdoses3 = 3 if maxdoses3==1 & n006_remark=="3针"
replace maxdoses3 = 3 if maxdoses3==1 & n006_remark=="三针"
replace maxdoses3 = 3 if maxdoses3==1 & n006_remark=="3"
replace maxdoses3 = 3 if maxdoses3==1 & n006_remark=="4针"
replace maxdoses3 = 3 if maxdoses3==1 & n006_remark=="应为三针"
replace maxdoses3 = 3 if maxdoses3==1 & n006_remark=="三针不包括加强针"
replace maxdoses3 = 3 if maxdoses3==1 & n006_remark=="三针，但是没有这个选项"
replace maxdoses3 = 3 if maxdoses3==1 & n006_remark=="三针，选项没这个"
replace maxdoses3 = 3 if maxdoses3==1 & usubstr(n006_remark,-8,8)=="不算加强针是三针"
replace maxdoses3 = 3 if maxdoses3==1 & usubstr(n006_remark,-8,8)=="不算加强针是3针"
replace maxdoses3 = 3 if maxdoses3==1 & usubstr(n006_remark,-8,8)=="不算加强针为3针"
replace maxdoses3 = 3 if maxdoses3==1 & usubstr(n006_remark,-8,8)=="不算加强针打3针"
replace maxdoses3 = 3 if maxdoses3==1 & usubstr(n006_remark,-7,7)=="不算加强针3针"
replace maxdoses3 = 3 if maxdoses3==1 & usubstr(n006_remark,-7,7)=="不算加强针三针" 
replace maxdoses3 = 3 if maxdoses3==1 & usubstr(n006_remark,-7,7)=="不算加强是三针" 
replace maxdoses3 = 3 if maxdoses3==1 & usubstr(n006_remark,-4,4)=="应为3针"
replace maxdoses3 = 3 if maxdoses3==1 & usubstr(n006_remark,-9,9)=="不包括加强针是3针"
replace maxdoses3 = 3 if maxdoses3==1 & usubstr(n006_remark,1,3)=="不包括"
replace maxdoses3 = 3 if maxdoses3==1 & ustrpos(n006_remark,"不算加强针打了三针")>0
replace maxdoses3 = 3 if maxdoses3==1 & ustrpos(n006_remark,"三针，没有加强针")>0
replace maxdoses3 = 3 if maxdoses3==1 & ustrpos(n006_remark,"总共打四针，前三针为基础针")>0
replace maxdoses3 = 3 if maxdoses3==1 & ustrpos(n006_remark,"不算加强针打了3针")>0
replace maxdoses3 = 3 if maxdoses3==1 & ustrpos(n006_remark,"第三针不是加强")>0
replace maxdoses3 = 3 if maxdoses3==1 & ustrpos(n006_remark,"3针都是基础针")>0
replace maxdoses3 = 3 if maxdoses3==1 & ustrpos(n006_remark,"3针基础，1针加强")>0
replace maxdoses3 = 3 if maxdoses3==1 & ustrpos(n006_remark,"安徽")>0
replace n006 = 3 if maxdoses3==3

keep primkey n004 n005 n006 
ren n004 vax_doses
recode vax_doses (0=0) (1/10=1), gen(vaxed)
replace vax_doses = 4 if inrange(vax_doses,4,100)
lab def vaxed 0 "unvaccinated" 1 "vaxed", replace
lab val vaxed vaxed
lab var vaxed "Vaccinated"
clonevar vaxed100 = vaxed
replace vaxed100 = vaxed100*100
lab def vax_doses  0 "0" 1 "1" 2 "2" 3 "3" 4 "4+", replace
lab val vax_doses vax_doses  
lab var vax_doses "Number of doses received"

ren n005 unvaxed_reason
ren n006 vax_doses_required
lab var vax_doses_required "Full vaccination requirement"

recode vax_doses (3/4=3), gen(vax_completion)
replace vax_completion = 2 if vax_completion==1 & vax_doses_required==1
replace vax_completion = 1 if inlist(vax_completion,1,2) & vax_doses_required==3
replace vax_completion = 3 if inlist(vax_completion,3) & vax_doses_required==3
lab def vax_completion  0 "unvaccinated" 1 "first vaccinated" 2 "fully vaccinated" 3 "booster vaccinated", replace
lab val vax_completion vax_completion
lab var vax_completion "Vaccination completion"
gen vax_completion1 = 100*inlist(vax_completion,0) if ~mi(vax_completion)
gen vax_completion2 = 100*inlist(vax_completion,1,2,3) if ~mi(vax_completion)
gen vax_completion3 = 100*inlist(vax_completion,2,3) if ~mi(vax_completion)
gen vax_completion4 = 100*inlist(vax_completion,3) if ~mi(vax_completion)

gen hhidpn = primkey, after(primkey)        // the ID corresponding to processed data

merge 1:1 hhidpn using $processed_data/demo_info_2022.dta
keep if inlist(_merge,2,3)
drop _merge 
replace primkey = hhidpn
replace ID = "id" + string(_n) if mi(ID)    // the ID in public released data

replace Age = Age +1 if indv_wave == 2021   // use age in 2022
recode Age (1/51=1)(52/54=2)(55/59=3)(60/64=4)(65/69=5)(70/74=6)(75/79=7)(80/200=8), gen(agegr8)
lab copy ag agegr8
lab def agegr8 1 "1-51" 2 "52-54", modify
lab val agegr8 agegr8
lab var agegr8 "Age group"
recode Age (1/51=1)(52/59=2)(60/69=3)(70/79=4)(80/200=5), gen(agegr5)
lab def agegr5 1 "1-51" 2 "52-59" 3 "60-69" 4 "70-79" 5 "80+"
lab val agegr5 agegr5
lab var agegr5 "Age group"
drop Age_Group  

lab var Edu_Group "Education group"
lab def edu_group 5 "High+", modify

recode Edu_Group (2/3=2) (4/5=3), gen(edugr3)
lab def edugr3 1 "Illiterate" 2 "Literate/Primary" 3 "Middle+"
lab val edugr3 edugr3
lab var Hukou "Hukou"

replace Female = 1 if Gender==2
replace Female = 0 if Gender==1


** ================================================================================
** 2. Data appending with weights, health-related variables, and vaccination in 2021
** ================================================================================

** ======== 2.1 Data appending ======================

// weights
merge 1:1 hhidpn using $cdata_22/COVID_Weight
drop if _merge==2
drop _merge

// ethnicity
merge 1:1 hhidpn using $processed_data/ethnicity.dta, keepus(minority)
drop if _merge==2
drop _merge

// health-related variables in 2018 and 2021/22
merge 1:1 ID using $varset_data/Health_Status/Var_Health_status_2018.dta, keepus(Chro_disease_any Cdisease_* ADL_needhelp_* IADL_needhelp_*  ADL_IADL_needhelp)
drop if _merge==2
drop _merge
drop Cdisease_11 Cdisease_12 //Emotional problems & Memory-related disease
ren (Cdisease_14 Cdisease_10) (Cdisease_10 Cdisease_14) // change position of Asthma and digestive disease
ren (Chro_disease_any Cdisease_* ADL_needhelp_* IADL_needhelp_* ADL_IADL_needhelp) w18_=
merge 1:1 hhidpn using $processed_data/Var_Health_status_2021.dta, keepus(Chro_disease_any Cdisease_* ADL_needhelp_* IADL_needhelp_* ADL_IADL_needhelp)
drop if _merge==2
drop _merge
merge 1:1 hhidpn using $processed_data/Var_Health_status_2022.dta, keepus(Chro_disease_any Cdisease_* ADL_needhelp_* IADL_needhelp_* ADL_IADL_needhelp) update replace
drop if _merge==2
drop _merge
drop Cdisease_11 Cdisease_12 Cdisease_13 Cdisease_14  //Emotional problems & Memory-related disease Parkinson disease Atrophy
ren (Cdisease_16 Cdisease_10) (Cdisease_10 Cdisease_16)
ren (Chro_disease_any Cdisease_* ADL_needhelp_* IADL_needhelp_* ADL_IADL_needhelp) w22_=

foreach w in 18 22 {
    recode w`w'_rating_health (997=.)
    replace w`w'_rating_health = w`w'_rating_health - 1
    drop w`w'_IADL_needhelp_phoning w`w'_ADL_needhelp_any w`w'_IADL_needhelp_any
    egen w`w'_ADL_IADL_items = rowtotal(w`w'_ADL_needhelp_* w`w'_IADL_needhelp_*), missing
    lab var w`w'_ADL_IADL_items "ADL/IADL dependencies"
    lab var w`w'_ADL_IADL_needhelp "ADL/IADL dependency"
    lab def w`w'_ADL_IADL_needhelp 0 "independent" 1 "dependent", replace
    lab val w`w'_ADL_IADL_needhelp ADL_IADL_needhelp
}

// vaccination data in the 2021 survey
merge 1:1 primkey using $fdata_21/Health_Insurance, keepus(eb001 eb035_1)
drop if _merge==2
drop _merge

replace eb001="1" if eb035_1=="新冠疫苗只接种了一针"  // a correction
charls_gen_mchoice eb001, items(1/4)    // convert answers with multiple selections to multiple binary variables
ren eb001_s1 vaxed_2021
drop eb001* eb035_1

gen vaxed_b21 = 0 if vaxed_2021==0 & vaxed==1
replace vaxed_b21 = 1 if vaxed_2021==1 & vaxed==1
gen vaxed_early = 0 if vaxed_2021==0 & vaxed==0 
replace vaxed_early = 1 if vaxed_2021==1 & vaxed==1
gen vaxed3 = 0 if vaxed==0 // balanced panel not required; vaxed_2021==0 & vaxed==0 before 2022/11/4
replace vaxed3 = 1 if vaxed_2021==0 & vaxed==1 
replace vaxed3 = 2 if vaxed_2021==1 & vaxed==1 
lab var vaxed3 "Vaxed (not yet/last year/earlier 0/1/2)"
lab def vaxed3 0 "still not vaccinated" 1 "last year" 2 "before last year", replace
lab val vaxed3 vaxed3
lab var vaxed3 "First vaccination time"

order w18_Cdisease_?, alpha
order w18_Cdisease_??, alpha after(w18_Cdisease_9)
order w22_Cdisease_?, alpha
order w22_Cdisease_??, alpha after(w22_Cdisease_9)

gen year = 2022

save $working_data/vax_sample_2022_raw, replace

** ======== 2.2 Sample selection (Figure 3) ======================

use $working_data/vax_sample_2022_raw, clear

    gen byte original = 0
    gl cond original
    gl cond $cond mi(INDV_weight_ad2)
    gl cond $cond mi(COVID_weight_ad2)
    gl cond $cond mi($weightvar)
    gl cond $cond Age<=51

    loc i = 4
    foreach c of global cond {
        drop if `c'
        // describing sample selection for Figure 3 in Methods
        count
        loc ++i
    }
    drop original

// region
recode provid (11 12 13 31 32 33 35 37 44 46 = 1) (14 34 36 41 42 43 = 2) (15 45 50/65 =3) (21 22 23 = 4), gen(region)
// 1=east, 2=central, 3=west, 4=northeast
lab var region "Region"  
lab def region 1 "East" 2 "Central" 3 "West" 4 "Northeast"
lab val region region

compress
save $working_data/vax_sample_2022, replace


** =================================================================
** 3. Description
** =================================================================

// Table 1: Study participants description

use $working_data/vax_sample_2022, clear
svyset [pw=$weightvar]

order w22_Cdisease_1 w22_Cdisease_7 w22_Cdisease_8 w22_Cdisease_3 w22_Cdisease_5 w22_Cdisease_4, before(w22_Cdisease_2)
order w18_Cdisease_1 w18_Cdisease_7 w18_Cdisease_8 w18_Cdisease_3 w18_Cdisease_5 w18_Cdisease_4, before(w18_Cdisease_2)

tab agegr5, matcell(freq)

foreach wgt in 1 $weightvar {        // the loop with (1) no weights (2) COVID-19 non-response adjusted weights
    foreach v of varlist Female Married vax_completion minority edugr3 Hukou region  ///
            w22_ADL_IADL_needhelp w22_Cdisease_* {
        tab `v' [aw=`wgt'], matcell(matall) 
        tab `v' agegr5 [aw=`wgt'], matcell(matage)
    }
}

// Figure 2: CHARLS statistics of vaccination rates


use $working_data/vax_sample_2022, clear 
svyset [pw=$weightvar]

    gen Age2021 = Age - 1 // adjust Age in 2021
    recode Age2021 (1/51=1)(52/59=2)(60/69=3)(70/79=4)(80/200=5), gen(agegr5_2021)
    lab def agegr5 1 "1-51" 2 "52-59" 3 "60-69" 4 "70-79" 5 "80+", replace
    lab val agegr5_2021 agegr5
    lab var agegr5 "Age group" 
    
    // first dose vaccination, 2021    
    recode vaxed3 (1=0) (2=100), gen(vaxed_21)
    svy: mean vaxed_21 if inrange(Age2021,52,.), over(agegr5_2021)

    // first dose vaccination, primary series completion, and booster completion by age group, 2022
    svy: mean vaxed100 vax_completion3 vax_completion4, over(agegr5)

// ODI Figure 1: A comparison with NHC's statistics

use $working_data/vax_sample_2022, clear
svyset [pw=$weightvar]

recode Age (60/69=1) (70/79=2) (80/200=3) (nonmiss=0), gen(agegr4)

    // CHARLS: vaccination, primary series completion, and booster completion by age group, 2022
    svy: mean vax_completion2-vax_completion4 if inrange(Age,60,.), over(agegr4) 



** =================================================================
** 4. Regression
** =================================================================


** ======== 4.1 variable preparation ======================

use $working_data/vax_sample_2022, clear 

** Checking the number of missing observations

    gl cond 
    gl cond $cond mi(Gender)
    gl cond $cond mi(Age)
    gl cond $cond mi(Edu_Group)
    gl cond $cond mi(Married)
    gl cond $cond mi(minority)
    gl cond $cond mi(Hukou)
    gl cond $cond mi(w22_Chro_disease_any)
    gl cond $cond mi(w22_ADL_IADL_needhelp)  

    loc i = 13
    foreach c of global cond {
        cou if `c'
        loc ++i
    }

// generate dependent variables
gen byte vaxed_fully = inlist(vax_completion,2,3)
gen byte vaxed_booster = inlist(vax_completion,3)

gen vaxed_fully1 = vaxed_fully
replace vaxed_fully1 = . if vaxed==0
gen vaxed_booster1 = vaxed_booster
replace vaxed_booster1 = . if vaxed_fully==0

gen vaxed_blced = vaxed if ~mi(vaxed3)
gen vaxed_fully_blced = vaxed_fully if ~mi(vaxed3)
gen vaxed_booster_blced = vaxed_booster if ~mi(vaxed3)

// a measure of comorbidity
egen w22_disease_number = rowtotal(w22_Cdisease_1-w22_Cdisease_10), missing
egen w18_disease_number = rowtotal(w18_Cdisease_1-w18_Cdisease_10), missing

// the household ID for standard error clustering
egen hhid1 = group(hhid)

save $working_data/vax_reg_2022, replace


** ======== 4.2 summary statistics ======================

use $working_data/vax_reg_2022, clear 

gl yvar vaxed vaxed_fully vaxed_booster
gl yvars1 vaxed_fully1 vaxed_booster1        // versus 1 dose or primary
gl yvars2 vaxed_blced vaxed_fully_blced vaxed_booster_blced        // versus 1 dose or primary
gl yvar1 vaxed3 
gl xvarsum1 Female Married minority i.agegr5 i.edugr3 i.Hukou i.region
gl xvarsum2_w22 w22_ADL_IADL_needhelp w22_Cdisease_* 
gl xvarsum2_w18 w18_ADL_IADL_needhelp w18_Cdisease_* 
gl xvarsum $xvarsum1 $xvarsum2_w22 $xvarsum2_w18
xi, noomit: tabstat $yvar $yvars1 $yvars2 $yvar1 $xvarsum, s(mean sd min max count) save c(stat)

// summary stats
mat sum = r(StatTotal)'

// imputation by creating a missing category
recode Edu_Group edugr3 Married minority w*Cdisease_* w*_disease_number w*ADL_IADL_needhelp w*ADL_IADL_items (missing = 99)
gen w22_disease_number_mi = mi(w22_disease_number)
gen w18_disease_number_mi = mi(w18_disease_number)

order w22_Cdisease_1 w22_Cdisease_7 w22_Cdisease_8 w22_Cdisease_3 w22_Cdisease_5 w22_Cdisease_4, before(w22_Cdisease_2)
order w18_Cdisease_1 w18_Cdisease_7 w18_Cdisease_8 w18_Cdisease_3 w18_Cdisease_5 w18_Cdisease_4, before(w18_Cdisease_2)
gl xvar1 i.Female i.Married i.minority i.Hukou i.agegr5 i.edugr3 i.region
gl xvar21_w22 i.w22_ADL_IADL_needhelp i.w22_Cdisease_*
gl xvar21_w18 i.w18_ADL_IADL_needhelp i.w18_Cdisease_*
gl xvar22_w22 i.w22_ADL_IADL_needhelp w22_disease_number i.w22_Cdisease_15 i.w22_Cdisease_16
gl xvar22_w18 i.w18_ADL_IADL_needhelp w18_disease_number i.w18_Cdisease_13 i.w18_Cdisease_14

** ======== 4.3 regressions ======================

est clear

// base outcome for mlogit
loc base_xvar21_w22_1 0
loc base_xvar21_w22_2 2
loc base_xvar22_w22_1 0
loc base_xvar22_w22_2 2
loc base_xvar21_w18_1 0
loc base_xvar21_w18_2 2
loc base_xvar22_w18_1 1
loc base_xvar22_w18_2 1

loc mlgtsd18    // no clustering due to a convergence issue
loc mlgtsd22 cluster hhid

foreach w in 22 18 {    // the loop with health-related variables either from 2021/22 or from 2018

// results in Table 2 (column 1/2/3): vaccination regression
// the results with health-related measures in 2018 are shown in ODI Table 2
    foreach y of global yvar {      // the loop with (1) vaccination (2) primary series completion (3) booster completion
        forv k = 1/2 {      // the loop with (1) binary indicators of chronic diseases (2) a comorbidity measure (ODI table 4)
            logit `y' $xvar1 ${xvar2`k'_w`w'}, vce(cluster hhid) or 
        }
    }

// results in Table 2 (column 4/5): change dependent variables (comparison groups) in logistic regression
    foreach y of global yvars1 {
        logit `y' $xvar1 ${xvar21_w22}, vce(cluster hhid) or 
    }

    
// results in Table 3: regression on vaccination timing
// the results with health-related measures in 2018 are shown in ODI Table 5
    forv k = 1/2 {          // the loop with (1) binary indicators of chronic diseases (2) a comorbidity measure
        loc baseo = `base_xvar2`k'_w`w'_1'
        mlogit $yvar1 $xvar1 ${xvar2`k'_w`w'}, vce(`mlgtsd`w'') rrr base(`baseo')

        // second-time regression to get confidence intervals for inverted RRR
        loc baseo = `base_xvar2`k'_w`w'_2'      
        mlogit $yvar1 $xvar1 ${xvar2`k'_w`w'}, vce(`mlgtsd`w'') rrr base(`baseo')
        
    }

// supplementary results (1): balanced panel in logistic regression
// ODI Table 3

    foreach y of global yvars2 {
        logit `y' $xvar1 ${xvar21_w22}, vce(cluster hhid) or 
    }


