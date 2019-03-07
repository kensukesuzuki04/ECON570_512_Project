clear
clc

delete MLE.txt
diary MLE.txt
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAXMIMIZING LOG LIKELIFOOD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Use HS
usehs = 0;

usegridn = 0;

% ISIC code
isic = 29;

% Loading Data and Parameters
A01_LoadData;

% lnk and lnn grid
A02_lnk_assignN_grid;

% First stage estimation
%B01_FirstStage_0221_onlyextensive_VFI;
B01_FirstStage

% Random grid on x
C02_RandomGrid;

% Productivity transition
C03_ProductivityTransition;

% Profit
C04_NonimporterProfit;

C05_ImporterProfit;


% %% Compute initial year probability
as_int = [-9 4 2];
% 
% x0=[as_int]';
% options = optimset('Display','iter','TolFun',1e-3,'TolX',1e-6, 'MaxIter',100000,'MaxFunEvals',100000);
% [para_as,fval,exitflag,output] = fminsearch(@F04loglf_firstyear,x0,options);
% 
% a0 = para_as(1)
% a1 = para_as(2)
% a2 = para_as(3)

%% Compute gamma distribution


gamma_int=[10 100]; % gammaIF, marginS

tic
disp(['******* estimation with fmin search / fminunc  ******* ']);

%S01startingvalue;
x0=[gamma_int, as_int];

options = optimset('Display','iter','TolFun',1e-4,'TolX',1e-8, 'MaxIter',100000,'MaxFunEvals',100000);
[fmin_parameter,fval,exitflag,output] = fminsearch(@F04_loglf,x0,options);
gammaIF = fmin_parameter(1);
gammaNS = fmin_parameter(2);
a0 = fmin_parameter(3);
a1 = fmin_parameter(4);
a2 = fmin_parameter(5);
par.gamma0=fmin_parameter(1:2);
gammaIF_fmin = fmin_parameter(1);
gammaNS_fmin = fmin_parameter(2);


% fmin unc
options = optimset('Display','iter','TolFun',1e-4,'TolX',1e-8, 'MaxIter',100000,'MaxFunEvals',100000);
[fminunc_parameter,fval,exitflag,output,grad,hessian] = fminunc(@F04_loglf,x0,options)
fminunc.std =  sqrt(diag(inv(hessian)))
fminunc.par =  fminunc_parameter

save R_MLEfminunc fminunc

fmintoc = toc




%% anneal

tic
for k=1:1
disp(['******* estimation with simulated annealing ******* ']);
%S01startingvalue;
x0=[fmin_parameter];

[anneal_parameter,fval_anneal] = anneal(@F04_loglf,x0); 

end
annealtoc = toc

save R_MLEpara fmin_parameter anneal_parameter

diary off


load R_MLEfminunc
load R_MLEpara

