% Main code to use MCMC runs

clear all
clc

% Use HS
usehs = 0;

% ISIC code
isic = 29;

% Loading Data and Parameters
A01_LoadData;

% lnk and lnn grid
A02_lnk_assignN_grid;

% First stage estimation
B01_FirstStage

% Random grid on x
C02_RandomGrid;

% Productivity transition
C03_ProductivityTransition;

% Profit
C04_NonimporterProfit;

C05_ImporterProfit;

% setting MH simulations
nrep = 25000; % number ofd repetitions
nparm = 5; % total number of parameters

% draw the random uniform noises for each iteration
um = log(rand(nrep,2));

% define the parameter vectors
par = zeros(nrep,nparm);
rate = zeros(nrep,1);

%% Setup the block of parameters and prior distribution parameters
% import fixed cost and sunk cost
niparm = 2;
% gammaIF, gammaNS, normal(0,1000)
mgammaFS = [0 0 ]; % mean
vgammasFS = [1000 1000]; % variance
priormi = [mgammaFS]; % prior mean
priorvi = [vgammasFS]; % prior var

% initial condition of importing
% a0 a1 a2, normal(0,100)
naparm=3;
priorma = [0 0 0 ]; % prior mean
priorva = [100 100 100]; % prior variance

% this is the end of the prior setting

%% steps for random work
stepi1 = ones(1,niparm/2)*.25; % for fixed cost
stepi2 = ones(1,niparm/2)*1; % for sunkcost
stepi = [stepi1 stepi2];
stepa = ones(1,naparm)*0.005;

%% Define initial starting values
istart = [3.1699 85.3917];
astart = [-6.1272 -0.6647 0.5056];
start = [istart astart];

%% Start the mcmc algorithm
par(1,:) = start;
rate(1) = 1;

parcurr = start; % current parameter

% keep acceptance rate
tacc = zeros(nrep,1); % ?
tacca = zeros(nrep,1); % for initial condition
tacci = zeros(nrep,1); % for fixed/sunk cost

record = zeros(nrep,4); % what this 4 means?

for i = 2:1:nrep
    %%% update paramerters %%%
    F05_updatepar; % construct [icurr acurr] from parcurr
    
    if i == 2
        % compute prior density
        priden =[normpdf(icurr,priormi,priorvi) normpdf(acurr,priorma,priorva)];
        % take log
        lpriden = log(priden);
        % log likelihood at start value
        ln= - F04_loglf(parcurr);
        % total 
        rllcurr = ln + sum(lpriden);
    end
    
    %%% update a %%%
    F05_updatepar;
    aprop = acurr + randn(1,naparm).*stepa; % proposal
    
    % prior
    priden =[normpdf(icurr,priormi,priorvi) normpdf(aprop,priorma,priorva)];
    % take log
    lpriden = log(priden);
    % update parprop
    parprop = [icurr aprop]; % change parameter for initial condition
    % log likelihood at start value
    ln= - F04_loglf(parprop);
    rllprop=ln + sum(lpriden);
    % update the parameter vector
    d = rllprop -rllcurr;
    % acceptance probability
    accprob = min(0,d); 
    
    if um(i,1)<=accprob; % accepted
        parcurr = parprop;
        tacca(i) = tacca(i-1) + 1;
        rllcurr = rllprop; % update
        display('proposal is accepted')
        display('initial condition a')
        [aprop]
    else
        parcurr = parcurr;
        tacca(i) = tacca(i-1);
    end
    
    %%% update fixed/sunk cost %%%
    F05_updatepar;
    iprop = icurr + randn(1,niparm).*stepi;
    %prior
    priden =[normpdf(iprop,priormi,priorvi) normpdf(acurr,priorma,priorva)];
    % take log
    lpriden = log(priden);
    % update parprop
    parprop = [iprop acurr]; % change parameter for initial condition
    % log likelihood at start value
    ln= - F04_loglf(parprop);
    rllprop=ln + sum(lpriden);
    % update the parameter vector
    d = rllprop -rllcurr;
    % acceptance probability
    accprob = min(0,d);
    
    
    if um(i,2)<=accprob; % accepted
        parcurr = parprop;
        tacci(i) = tacci(i-1) + 1;
        rllcurr = rllprop; % update
        display('proposal is accepted')
        display('import fixed/sunk cost')
        [iprop]
    else
        parcurr = parcurr;
        tacci(i) = tacci(i-1);
    end
    
    % save optimal for this iteration
    par(i,:) = parcurr;
    
    tacc(i) = min([tacca(i) tacci(i)]);
    rate(i) = tacc(i)/(i-1);
    disp('number of iteration')
    i
    disp('accept rate')
    rate(i)
    
    record(i,1) = rllcurr;
    record(i,2) = rllprop;
    record(i,3) = accprob;
    record(1,3) = um(i,2);
    
    save MCMCresult par rate record
end

    
    
