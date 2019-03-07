%run A01_LoadData

global sigmahat alpha0 alpha1 alpha2 alpha3 alpha4 ...
    betak betan omega rev_intcpt sigmaxi

%% STEP 1: Recovering Elasticty
% Reg eq: tvc = zeta * income_sales + nu

Elast.zetahat = (income_sales' * income_sales)\income_sales'*tvc;   % coefficient
Elast.resid = tvc - income_sales*Elast.zetahat;                   % residual
Elast.nuhat = (Elast.resid'*Elast.resid) / (N-1);                 % variance of error
Elast.var = Elast.nuhat * inv(income_sales' * income_sales);      % variance-covariance matrix for coefficients
Elast.tstat = Elast.zetahat / Elast.var^(1/2);                    % T-statistics (betanot = 0)
Elast.p = 1 - tcdf(Elast.tstat,N-1);                              % p-value
Elast.Rsq = 1- (Elast.resid'*Elast.resid) ...
            / (((tvc - mean(tvc)))'*((tvc - mean(tvc)))); % R-squared:1- SSE / SSTO (total sum of squared)
        
disp('----------------------------')
disp('STEP1 Estimaton')
Elast
disp('Elasticity of substitution: sigma = 1/(1-sigma)')
Para.sigmahat = 1/(1-Elast.zetahat)
disp(' ')

%% Step 2: Nonparametric Revenue Production Funcion Estimation
% Time dummy (2000 as reference)
D.T2001 = zeros(N,1);               % 2001 dummy
D.T2001(find(year==2001)) = 1;
D.T2002 = zeros(N,1);               % 2002 dummy
D.T2002(find(year==2002)) = 1;
D.T2003 = zeros(N,1);               % 2003 dummy
D.T2003(find(year==2003)) = 1;
D.T2004 = zeros(N,1);               % 2004 dummy
D.T2004(find(year==2004)) = 1;
D.T2005 = zeros(N,1);               % 2005 dummy
D.T2005(find(year==2005)) = 1;
D.T2006 = zeros(N,1);               % 2006 dummy
D.T2006(find(year==2006)) = 1;
D.T = [D.T2001, D.T2002, D.T2003, D.T2004, D.T2005, D.T2006];

clear  lnn2 lnn3 kn kn2 k2n nx nx2 n2x

grid_lnn2_ = grid_lnn1_ .^2;
grid_lnn3_ = grid_lnn1_ .^3;
kn = lnk1 .* grid_lnn1_;
kn2 = lnk1 .* grid_lnn2_;
k2n = lnk2 .* grid_lnn1_;
nx = grid_lnn1_ .* lnx1;
nx2 = grid_lnn1_ .* lnx2;
n2x = grid_lnn2_ .* lnx1;

grid_lnnhs2_ = grid_lnnhs1_ .^2;
grid_lnnhs3_ = grid_lnnhs1_ .^3;
knhs = lnk1 .* grid_lnnhs1_;
knhs2 = lnk1 .* grid_lnnhs2_;
k2nhs = lnk2 .* grid_lnnhs1_;
nhsx = grid_lnnhs1_ .* lnx1;
nhsx2 = grid_lnnhs1_ .* lnx2;
nhs2x = grid_lnnhs2_ .* lnx1;


if usegridn == 0
    grid_lnn1_ = lnn1;
    grid_lnn2_ = lnn1 .^2;
    grid_lnn3_ = lnn1 .^3;
    kn = lnk1 .* lnn1;
    kn2 = lnk1 .* (lnn1.^2);
    k2n = lnk2 .* lnn1;
    nx = lnn1 .* lnx1;
    nx2 =lnn1 .* lnx2;
    n2x =(lnn1.^2) .* lnx1;
else
end

% Stat.X = [ones(N,1), D.T, lnk1, lnk2, lnk3, lnx1, lnx2, lnx3, grid_lnnhs1_, grid_lnnhs2_ grid_lnnhs3_,...
%             kx, kx2, k2x, knhs, knhs2, k2nhs, nhsx, nhsx2, nhs2x];

Stat.X = [ones(N,1), D.T, lnk1, lnk2, lnk3, lnx1, lnx2, lnx3, grid_lnn1_, grid_lnn2_ grid_lnn3_,...
            kx, kx2, k2x, kn, kn2, k2n, nx, nx2, n2x];

Stat.phihat = (Stat.X'*Stat.X)\ (Stat.X'*lnrev); % estimated coefficient
Stat.resid = lnrev - Stat.X* Stat.phihat;        % residual
Stat.nuhat = (Stat.resid'*Stat.resid) / (N-size(Stat.X,2));   % variance of error (unbiased estimate for variance)
Stat.var = Stat.nuhat * inv(Stat.X' * Stat.X);   % variance-covariance matrix for coefficients
Stat.std = diag(Stat.var).^(1/2);
Stat.tstat = Stat.phihat ./ Stat.std;       % T-statistics (betanot = 0)
Stat.p = tcdf(-abs(Stat.tstat),N-size(Stat.X,2))  + 1 -  tcdf(abs(Stat.tstat),N-size(Stat.X,2));
Stat.Rsq = 1- (Stat.resid'*Stat.resid) ...
            / (((lnrev - mean(lnrev)))'*((lnrev - mean(lnrev)))); % R-squared:1- SSE / SSTO (total sum of squared)

Stat.phiselect = Stat.phihat([1 8:length(Stat.phihat)]);
Stat.pselect   = Stat.p([1 8:length(Stat.phihat)]);
Stat.stdselect   = Stat.std([1 8:length(Stat.phihat)]);
Stat.tstatselect   = Stat.tstat([1 8:length(Stat.phihat)]);
Stat.indexselect = {'Const'; 'lnk1'; 'lnk2'; 'lnk3'; 'lnx'; 'lnx2'; 'lnx3'; ...
                    'lnn1';  'lnn2'; 'lnn3'; 'kx'; 'kx2'; 'k2x'; ...
                    'kn';  'kn2'; 'k2n'; 'nx'; 'nx2'; 'n2x'; };
% Stat.indexselect = {'Const'; 'lnk'; 'lnk2'; 'lnk3'; 'lnx'; 'lnx2'; 'lnx3'; ...
%                     'n';  'n2'; 'n3'};
table(Stat.indexselect, Stat.phiselect, Stat.tstatselect, Stat.pselect)

disp('----------------------------')
disp('STEP2: Static Estimaton')
Stat
disp(' ')


%% Step 3: NLS for Evoluition for productivity
Evol.varphi = Stat.phihat(8:length(Stat.phihat),1); % coefficeints of nonparametric function h(varphi)
Evol.data = Stat.X(:,8:(length(Stat.phihat)));
Evol.phihat = Evol.data * Evol.varphi;

% generate lagged phi
Evol.lagphihat = zeros(T*Nf,1);
Evol.lagphihat(2:T*Nf) = Evol.phihat(1:(T*Nf - 1));
Evol.lagphihat = Evol.lagphihat.*notfirst;

% generate lagged k
Evol.laglnk = zeros(T*Nf,1);
Evol.laglnk(2:T*Nf) = lnk1(1:(T*Nf - 1));
Evol.laglnk = Evol.laglnk.*notfirst;

% generate lagged relative importance of import
Evol.laglnn_ = zeros(T*Nf,1);
Evol.laglnn_(2:T*Nf) = grid_lnn1_(1:(T*Nf - 1));
Evol.laglnn_ = Evol.laglnn_.*notfirst;

% generate lagged relative importance of import
Evol.laglnnhs_ = zeros(T*Nf,1);
Evol.laglnnhs_(2:T*Nf) = grid_lnnhs1_(1:(T*Nf - 1));
Evol.laglnnhs_ = Evol.laglnnhs_.*notfirst;

% generate lagged import dummy
Evol.lagimp_dummy = zeros(T*Nf,1);
Evol.lagimp_dummy(2:T*Nf) = imp_dummy(1:(T*Nf - 1));
Evol.lagimp_dummy = Evol.lagimp_dummy.*notfirst;

% NLS
Evol.alpha = ones(5,1);     % alpha0, alpha1, ..., alpha4
Evol.beta = ones(2,1);       % etak etan
Evol.para = [Evol.alpha;Evol.beta];

Evol.para_int = ones(length(Evol.para),1);

IndexwoF = find(notfirst==1);   % index without first year
NwoF = length(IndexwoF);        % Number of observation

% 
% NLSevol = @(para) F01_NLSevol (para, ...  % parameter
%     Evol.lagphihat(IndexwoF), ...       % lagged phihat
%     Evol.laglnk(IndexwoF), ...          % Lagged lnk this is continuous
%     Evol.laglnnhs_(IndexwoF), ...       % Lagged lnn
%     Evol.lagimp_dummy(IndexwoF), ...    & lagged import dummy
%     lnk1(IndexwoF), ...                 % lnk this is continuous
%     grid_lnnhs1_(IndexwoF),...          % lnn this is exogenous (use bin)
%     Evol.phihat(IndexwoF), ...          % phihat
%     Para.sigmahat, ...                  % elasticity
%     NwoF);                              % Number of obs (wo first year)

NLSevol = @(para) F01_NLSevol (para, ...  % parameter
    Evol.lagphihat(IndexwoF), ...       % lagged phihat
    Evol.laglnk(IndexwoF), ...          % Lagged lnk
    Evol.laglnn_(IndexwoF), ...         % Lagged lnn: this is exogenous (use bin)
    Evol.lagimp_dummy(IndexwoF), ...    & lagged import dummy
    lnk1(IndexwoF), ...                 % lnk
    grid_lnn1_(IndexwoF),...            % lnn: this is exogenous (use bin)
    Evol.phihat(IndexwoF), ...          % phihat
    Para.sigmahat, ...                  % elasticity
    NwoF);                              % Number of obs (wo first year)
% anonumous function of parameter (used for lsqnonlin)

options = optimoptions('lsqnonlin','Display','iter');
[Evol.para,Evol.SSR,Evol.resids,Evol.exitflag,Evol.output,Evol.lambda,Evol.J]...
    = lsqnonlin(NLSevol, Evol.para_int, [], [], options);
% para: solution to parameter
% SSR: sum of squared residual (scalar)
% res: residual (N x 1 vector)
% exitflag: excit condition
% output: information about optimization
% lambda: Langange multipliers at the solution (para)
% J: Jacobian of the function at the solution

Evol.xis = Evol.SSR / (NwoF - size(Evol.para,1));   % variance of error (unbiased estimate for variance)
Evol.cov = inv(Evol.J'*Evol.J)*var(Evol.resids);       % variance covariance matrix of estimates
Evol.std = diag(Evol.cov).^(1/2);                   % standard error
Evol.tstat = Evol.para ./ Evol.std;                 
Evol.p = tcdf(-abs(Evol.tstat),N-size(Evol.para,1))  + 1 -  tcdf(abs(Evol.tstat),N-size(Evol.para,1));

Evol.Rsq = 1- (Evol.resids'*Evol.resids) ...
            / (((Evol.phihat(IndexwoF) - mean(Evol.phihat(IndexwoF))))'*((Evol.phihat(IndexwoF) - mean(Evol.phihat(IndexwoF))))); % R-squared:1- SSE / SSTO (total sum of squared)
Evol.resid =Evol.resids  / (1-Para.sigmahat);
Para.sigmaxi = var(Evol.resid)^(1/2)
        
Evol.paraindex = {'alpha0';'alpha1';'alpha2';'alpha3';'alpha4'; 'beta_k';'beta_n'};
table(Evol.paraindex, Evol.para, Evol.std, Evol.tstat, Evol.p)

% Global variables
sigmahat = Para.sigmahat;
alpha0 = Evol.para(1);
alpha1 = Evol.para(2);
alpha2 = Evol.para(3);
alpha3 = Evol.para(4);
alpha4 = Evol.para(5);
betak =  Evol.para(6);
betan =  Evol.para(7);
sigmaxi = Para.sigmaxi;

%% Estimated Productivity

%omega = (-1)/(1-Para.sigmahat)*Evol.phihat + Evol.para(1,1)*lnk1 + Evol.para(2,1)*grid_lnn1_;
omega = (-1)/(1-Para.sigmahat)*Evol.phihat + betak*lnk1 + betan*grid_lnn1_;

%% Revenue intercept

rev_intcpt = Stat.phihat(1,1) + mean(Stat.phihat(2:7));

%% Export to CSV
clear exp
R_EstimatedOmega = table(id,year,imp_dummy,omega,lnk1, grid_lnk1);
writetable(R_EstimatedOmega)

save R_FirstStage Evol Stat 