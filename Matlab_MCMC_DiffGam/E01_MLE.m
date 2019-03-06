clear
clc

delete MLE.txt
diary MLE.txt
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAXMIMIZING LOG LIKELIFOOD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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

%READ IN ESTIMATED/STARTING VALUES OF PARAMETERS
SchEst=zeros(1,2);
SchLlh=zeros(1,1);
save SchEst SchEst
save SchLlh SchLlh

%% Compute gamma distribution

as_int = [-9 4 2];
gamma_int=[10*ones(1,numkgrid) 100*ones(1,numkgrid)]; % gammaIF, marginS

tic
for k=1:1
disp(['******* estimation ' num2str(k,'%10.0f') ' ******* ']);
%S01startingvalue;
x0=[gamma_int as_int];

options = optimset('Display','iter','TolFun',1e-4,'TolX',1e-8, 'MaxIter',100000,'MaxFunEvals',100000);
[parameter,fval,exitflag,output] = fminsearch(@F04loglf,x0,options);

parameter=parameter'
fval
par.gamma0=parameter(1:2);

gammaIF_fmin = parameter(1,1);
gammaNS_fmin = parameter(1,2);

save SchEst SchEst
save SchLlh SchLlh

end
fmintoc = toc

gammaIF = parameter(1);
gammaNS = parameter(2);
a0 = parameter(3);
a1 = parameter(4);
a2 = parameter(5);

tic
for k=1:1
disp(['******* estimation with simulated annealing ******* ']);
%S01startingvalue;
x0=[parameter]';

[parameter_anneal,fval_anneal] = anneal(@F04loglf,x0); 

parameter_anneal=parameter_anneal'
fval_anneal

end
annealtoc = toc

diary off