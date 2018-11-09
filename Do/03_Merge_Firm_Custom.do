clear all

use "C:\Users\KensukeSuzuki\Documents\2018\Firm level data\Chinese data\Cleaned\import_aggreg_allyears.dta"

drop if plant_id == ""
drop if year == .
save "C:\Users\KensukeSuzuki\Documents\2018\Firm level data\Chinese data\Cleaned\import_aggreg_allyears_conc.dta", replace

clear

use "C:\Users\KensukeSuzuki\Documents\2018\Firm level data\Chinese data\Cleaned\firm_allyears.dta", replace

sort plant_id year
by plant_id year : gen dup = cond(_N==1,0,_n)
tab dup

drop if dup> 0

merge 1:1 plant_id year using "C:\Users\KensukeSuzuki\Documents\2018\Firm level data\Chinese data\Cleaned\import_aggreg_allyears_conc.dta"

drop num_plant_id
order year plant_id party_id

drop if year < 2000

save "C:\Users\KensukeSuzuki\Documents\2018\Firm level data\Chinese data\Cleaned\firm_imp_comb.dta", replace
