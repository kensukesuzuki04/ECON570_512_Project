clear all

// Change directory
cd "C:\Users\KensukeSuzuki\Box Sync\2018\03 Fall\ECON570 DevEcon\Minipaper\ECON570_512_Project\Do"

local year = 1996
use "C:\Users\KensukeSuzuki\Documents\2018\Firm level data\Chinese data\Firms\unbalanced_firm_`year'.dta"



forval year = 1996(1)2006{
use "C:\Users\KensukeSuzuki\Documents\2018\Firm level data\Chinese data\Firms\unbalanced_firm_`year'.dta"

keep year firm_code industry_code isic_rv3 ownership_id start_year ///
	paidin_capital state_capital col_capital legal_person_capital ///
	private_capital hkt_capital for_capital inventory finished_p_inventory ///
	long_invest sales_income sales_profit t_profit ///
	inter_input export employee value_added r_d_cost isic4_update ///
	start_year

destring firm_code, generate(num_firm_code) force /*generate a new variable*/
rename firm_code plant_id
rename num_firm_code num_plant_id /* to match with custom data*/

destring industry_code, generate(num_industry_code) force /*generate a new variable*/
destring isic_rv3, generate(num_isic_rv3) force /*generate a new variable*/
destring isic4_update, generate(num_isic4_update) force /*generate a new variable*/

save "C:\Users\KensukeSuzuki\Documents\2018\Firm level data\Chinese data\Cleaned\firm_`year'.dta", replace
clear
}

// 1 SOE 
// 2 Collective and others
// 3 Foreign
// 4 HMT: Hong Kong Macao Taiwan?
// 5 private 

// Append
use "C:\Users\KensukeSuzuki\Documents\2018\Firm level data\Chinese data\Cleaned\firm_1996.dta"
forval year = 1997(1)2006{
append using "C:\Users\KensukeSuzuki\Documents\2018\Firm level data\Chinese data\Cleaned\firm_`year'.dta"
}
save "C:\Users\KensukeSuzuki\Documents\2018\Firm level data\Chinese data\Cleaned\firm_allyears.dta", replace
