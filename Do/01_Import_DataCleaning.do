clear all

// Change directory
cd "C:\Users\KensukeSuzuki\Box Sync\2018\03 Fall\ECON570 DevEcon\Minipaper\ECON570_512_Project\Do"


forval year = 2000(1)2006{
foreach month in "01" "02" "03" "04" "05" "06" "07" "08" "09" 10 11 12 { 

use "C:\Users\KensukeSuzuki\Documents\2018\Firm level data\Chinese data\Imports\Raw_import_`year'`month'.dta"

order year shipment_month party_id hs_id hs_id_num origin_id value quantity unit_id ownership_id ex_rate
keep year shipment_month party_id hs_id hs_id_num origin_id value quantity unit_id ownership_id ex_rate

// ownership id
// 2006
// 1: State owned
// 2: 
// 3: Foreign
// 4: Joint-venture
// 5: Private-domestic
// 2000
// 1: State owned
// 2: Collective owned
// 3: Foreign
// 4: Joint-venture
// 5: Private-domestic

destring hs_id, generate(num_hs_id) force /*generate a new variable*/
list hs_id if num_hs_id>=.	/*will show which have nonnumeric*/
drop if num_hs_id>=.
drop hs_id hs_id_num
rename num_hs_id hs_code_

destring party_id, generate(num_plant_id) force /*generate a new variable*/
list party_id if num_plant_id>=.	/*will show which have nonnumeric*/
//drop if num_plant_id>=.
//drop party_id
//rename num_party_id plant_id
rename party_id plant_id
format num_plant_id %10.0f

rename value value_
rename quantity quantity_
rename origin_id origin_id_

egen valuetemp = sum(value_), by(plant_id hs_code_ origin_id)
egen quantitytemp = sum(quantity_), by(plant_id hs_code_ origin_id)

bysort plant_id hs_code_ origin_id: gen keep = _n
sort hs_code_ origin_id, stable

keep if keep == 1
drop keep
drop quantity_ value_ 
rename quantitytemp quantity_
rename valuetemp value_

save "C:\Users\KensukeSuzuki\Documents\2018\Firm level data\Chinese data\Cleaned\import_hs_orig_`year'`month'.dta", replace

// Summing over country of origin
egen quantitytemp = sum(quantity), by(plant_id hs_code_)
egen valuetemp = sum(value_), by(plant_id hs_code_)
bysort plant_id hs_code_: gen keep = _n

// Compute the number of source countries
egen maxkeep =  max(keep), by(plant_id hs_code_)
keep if maxkeep == keep
drop keep origin_id
drop quantity_ value_ 

rename maxkeep num_origin_

rename quantitytemp quantity_
rename valuetemp value_

save "C:\Users\KensukeSuzuki\Documents\2018\Firm level data\Chinese data\Cleaned\import_hs_`year'`month'.dta", replace


// Summing over hs codes
egen valuetemp = sum(value_), by(plant_id)
bysort plant_id : gen keep = _n

// Compute the number of imported inputs
egen maxkeep =  max(keep), by(plant_id)
keep if maxkeep == keep
drop keep
drop  value_ 

rename maxkeep num_import_
rename valuetemp value_

drop hs_code 

save  "C:\Users\KensukeSuzuki\Documents\2018\Firm level data\Chinese data\Cleaned\import_`year'`month'.dta", replace
clear
}
}
//

// Combine monthly data
clear 
forval year = 2000(1)2006{

use  "C:\Users\KensukeSuzuki\Documents\2018\Firm level data\Chinese data\Cleaned\import_`year'01.dta"

foreach month in "02" "03" "04" "05" "06" "07" "08" "09" 10 11 12 { 
append using "C:\Users\KensukeSuzuki\Documents\2018\Firm level data\Chinese data\Cleaned\import_`year'`month'.dta"
}

drop quantity_

// Summing over months
egen valuetemp = sum(value_), by(plant_id)
bysort plant_id: gen keep = _n
keep if keep == 1
drop keep shipment_month
drop value_ 
rename valuetemp value_

rename plant_id party_id
merge 1:1 party_id using "C:\Users\KensukeSuzuki\Documents\2018\Firm level data\Chinese data\Custom_Firm_Concordance\ID_concordance`year'_adj.dta"
rename _merge merge_conc
rename firm_code plant_id

order year party_id plant_id num_plant_id

save  "C:\Users\KensukeSuzuki\Documents\2018\Firm level data\Chinese data\Cleaned\import_aggreg_`year'.dta", replace


clear
}
//


// Combine yearly data
clear
use  "C:\Users\KensukeSuzuki\Documents\2018\Firm level data\Chinese data\Cleaned\import_aggreg_2000.dta"

forval year = 2001(1)2006{
append using "C:\Users\KensukeSuzuki\Documents\2018\Firm level data\Chinese data\Cleaned\import_aggreg_`year'.dta"
}

save  "C:\Users\KensukeSuzuki\Documents\2018\Firm level data\Chinese data\Cleaned\import_aggreg_allyears.dta", replace


/////


