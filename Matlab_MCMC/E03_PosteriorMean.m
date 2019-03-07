% Get posterior mean

% specify burn-in period
burnin = 1000;

% Load result --- posterior draws
load R_MCMCresult

% Drop bunin period
st_par = par(burnin+1:length(par),:);

% posterior mean
ps_mean = mean(st_par)
ps_std = std(st_par)

