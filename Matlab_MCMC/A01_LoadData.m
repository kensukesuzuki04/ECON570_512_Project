%Load the data constructed from STATA
clear all
clc

%Load data 
data= readtable('C:\Users\KensukeSuzuki\Box Sync\2018\03 Fall\ECON570 DevEcon\Minipaper\ECON570_512_Project\isic29_priv_balanced.csv');
%data = readtable('C:\Users\KensukeSuzuki\Documents\2018\Firm level data\Chinese data\Cleaned\isic29_priv_balanced.csv');

global id year tvc income_sales N lnrev Nhs lnnhs1 lnk1 lnx1 lnn1 lagimp_dummy ...
       T Nf notfirst imp_dummy ...
       numkgrid numprodgrid scale ...
       delta ...
       intN S unique_id lastyear
S = 1;
   
% data
id = data.id;
year = data.year;
tvc = data.tvc;
income_sales = data.income_sales;
lnrev = data.log_rev;
lnk1 = data.lnk1;
lnk2 = data.lnk2;
lnk3 = data.lnk3;
lnx1 = data.lnx1;
lnx2 = data.lnx2;
lnx3 = data.lnx3;
lnn1 = data.lnn1;
lnn2 = data.lnn2;
lnn3 = data.lnn3;
kx  = data.crosskx;
kx2 = data.crosskx2;
k2x = data.crossk2x;
kn  = data.crosskn;
kn2 = data.crosskn2;
k2n = data.crossk2n;
nx  = data.crossnx;
nx2 = data.crossnx2; 
n2x = data.crossn2x;
T = 7;
Nf = length(id)/T;
notfirst = ones(T*Nf,1);
notfirst(find(year==2000),1) = 0;
lastyear = zeros(T*Nf,1);
lastyear(find(year==2006),1) = 1;
notfirstsecond = ones(T*Nf,1);
notfirstsecond(find(year==2000 | year == 2001),1) = 0;
N = length(id);
imp_dummy = data.import_dummy;

Nhs = data.num_hs;
Nhs(isnan(Nhs)) = 0;
Nhs = Nhs + 1;
lnnhs1 = log(Nhs);

unique_id = unique(id);

% Number of grid
numkgrid = 8;
numprodgrid = 100;

% Scaling
scale = 1e-3;

% Discount rate
delta = 0.8;

%# of quadrature points for numerical integration
intN = 30;

% generate lagged import dummy
lagimp_dummy = zeros(T*Nf,1);
lagimp_dummy(2:T*Nf) = imp_dummy(1:(T*Nf - 1));
lagimp_dummy = lagimp_dummy.*notfirst;
lagimp_dummy(find(notfirst==0)) = NaN;
