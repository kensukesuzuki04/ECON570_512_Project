clear all

// Change directory
cd "C:\Users\KensukeSuzuki\Box Sync\2018\03 Fall\ECON570 DevEcon\Minipaper\ECON570_512_Project\Do"


forval year = 2000(1)2006{
foreach month in "01" "02" "03" "04" "05" "06" "07" "08" "09" 10 11 12 { 

use "C:\Users\KensukeSuzuki\Documents\2018\Firm level data\Chinese data\Imports\Raw_import_`year'`month'.dta"

order year shipment_month party_id hs_id hs_id_num origin_id value quantity unit_id ownership_id ex_rate
keep year shipment_month party_id hs_id hs_id_num origin_id value quantity unit_id ownership_id ex_rate

destring hs_id, generate(num_hs_id) force /*generate a new variable*/
list hs_id if num_hs_id>=.	/*will show which have nonnumeric*/
drop if num_hs_id>=.
drop hs_id hs_id_num
rename num_hs_id hs_code_

destring party_id, generate(num_party_id) force /*generate a new variable*/
list party_id if num_party_id>=.	/*will show which have nonnumeric*/
drop if num_party_id>=.
drop party_id
rename num_party_id plant_id
format plant_id %10.0f

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

}
}
//

forval year = 2000(1)2006{

use  "C:\Users\KensukeSuzuki\Documents\2018\Firm level data\Chinese data\Cleaned\import_`year'01.dta"

foreach month in "02" "03" "04" "05" "06" "07" "08" "09" 10 11 12 { 
append using "C:\Users\KensukeSuzuki\Documents\2018\Firm level data\Chinese data\Cleaned\import_aggregate_select_`year'`month'.dta"
}



save  "C:\Users\KensukeSuzuki\Documents\2018\Firm level data\Chinese data\Cleaned\import_aggregate_select_`year'.dta", replace

}

