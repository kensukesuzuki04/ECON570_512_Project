%Load the data constructed from STATA

BaseName = 'C:\Users\KensukeSuzuki\Documents\2018\Firm level data\Chinese data\Cleaned\priv_balanced_isic';
FileName = [BaseName,num2str(isic),'.csv'];

%Load data 
data = readtable(FileName);

global id year tvc income_sales N lnrev Nhs lnnhs1 lnk1 lnx1 lnn1 lagimp_dummy ...
       T Nf notfirst imp_dummy ...
       numkgrid numprodgrid scale ...
       delta ...
       intN S unique_id lastyear usehs isic
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

Nhs = data.num_hs;
Nhs(isnan(Nhs)) = 0;
Nhs = Nhs + 1;
lnnhs1 = log(Nhs);

kx  = data.crosskx;
kx2 = data.crosskx2;
k2x = data.crossk2x;

if usehs == 0
    lnn1 = data.lnn1;
    lnn2 = data.lnn2;
    lnn3 = data.lnn3;
    
    kn  = data.crosskn;
    kn2 = data.crosskn2;
    k2n = data.crossk2n;
    nx  = data.crossnx;
    nx2 = data.crossnx2;
    n2x = data.crossn2x;
else
    lnn1 = lnnhs1 ;
    lnn2 = lnnhs1.^2;
    lnn3 = lnnhs1.^2;

    kn  = lnk1 .* lnn1;
    kn2 = lnk1 .* lnn2;
    k2n = lnk2 .* lnn1;
    nx  = lnn1 .* lnx1;
    nx2 = lnn1 .* lnx2;
    n2x = lnn1 .* lnx1;
end



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
