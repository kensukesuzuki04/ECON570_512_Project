% Monte Carlo Simulation

clear all

global Nf S T scale  notfirst numprodgrid ...
    xgrid lnk1 grid_lnk1 grid_lnn1 grid_lnnhs1 grid_lnn1_cand ...
    grid_kindex  omega ...
    piN piI TxN TxI ...
    imp_dummy lagimp_dummy ...
    numkgrid ...
    sigmaxi sigmahat delta ...
    betak betan alpha0 alpha1 alpha2 alpha3 alpha4 ...
    TxIobs TxNobs ...
    rev_intcpt ...
    a0 a1 a2 a3

load R_MCMCresult
E03_PosteriorMean
Sim.T = 20;

useMCMC = 1;
usehs = 0;
usegridn = 1;
isic = 29;

if useMCMC == 1
    Sim.gammaIF = ps_mean(1:8);
    Sim.gammaNS = ps_mean(9:16);
else
    Sim.gammaIF = fmiunc.par(1);
    Sim.gammaNS = fminunc.par(2);
end


%%
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


%% VFI

gammaIF = Sim.gammaIF;
gammaNS = Sim.gammaNS;
D01_VFI


%%

Sim.year = repmat([1:1:Sim.T]',Nf,1);
Sim.k = zeros(Nf*Sim.T, 1)
Sim.kind = zeros(Nf*Sim.T, 1)
Sim.n_cand = zeros(Nf*Sim.T, 1)
Sim.omega = zeros(Nf*Sim.T, 1)
Sim.imp = 999*ones(Nf*Sim.T, 1)
Sim.grid_kindex = zeros(Nf*Sim.T, 1)
Sim.FYgrid_kindex = zeros(Nf,1)

% assign capital (grid) and candidate n (grid) from data
for i = 1:Nf
    Sim.k((i-1)*Sim.T+1:i*Sim.T) = grid_lnk1((i-1)*T+1,1);
    Sim.grid_kindex((i-1)*Sim.T+1:i*Sim.T) = grid_kindex((i-1)*T+1,1);
    Sim.FYgrid_kindex = Sim.grid_kindex(find(Sim.year==1));
    Sim.n_cand((i-1)*Sim.T+1:i*Sim.T) = grid_lnn1_cand((i-1)*T+1,1); % grid_lnn1_cand assumes constant capital stock
end

% assing first year productivity
Sim.omega(find(Sim.year==1)) = omega(find(year==2000));

% use first year import status
Sim.imp(find(Sim.year==1)) = imp_dummy(find(year==2000));


% compute next period productivity (based on today's productivity and import status)
for t = 2:Sim.T
    
    seed = 19900904*t;
    rng(seed)
    % period t productiviy
    shock = normrnd(0,sigmaxi,[Nf,1]);
    Sim.omega(find(Sim.year==t)) = alpha0 + ...
        + alpha1*Sim.omega(find(Sim.year==(t-1))) ...
        + alpha2*Sim.omega(find(Sim.year==(t-1)) ).^2....
        + alpha3*Sim.omega(find(Sim.year==(t-1)) ).^3....
        + alpha4*Sim.imp(find(Sim.year==(t-1))) ...
        + shock;
    
    Sim.cost = zeros(Nf, 1);
    
    for cap = 1:numkgrid
       
    % draw relevant cost (depending on t-1 import status)
    gammaIFcap = Sim.gammaIF(cap);
    gammaNScap = Sim.gammaNS(cap);
    Sim.cost(find(Sim.FYgrid_kindex==cap)) =  Sim.imp(find(Sim.year==(t-1) & Sim.grid_kindex==cap)) .* exprnd(gammaIFcap,[length(Sim.imp(find(Sim.year==(t-1) & Sim.grid_kindex==cap))),1]) ... % case of importer
        + (1-Sim.imp(find(Sim.year==(t-1) & Sim.grid_kindex==cap))) .* exprnd(gammaNScap,[length(Sim.imp(find(Sim.year==(t-1) & Sim.grid_kindex==cap))),1]); % case of  non importer
    end
    
    % Profit at t
    SimpiI = exp(rev_intcpt +(1-sigmahat)*(betak*Sim.k(find(Sim.year==(t))) + betan * Sim.n_cand(find(Sim.year==(t))) - Sim.omega(find(Sim.year==(t))) )) ...
        *(1/sigmahat)*scale; % (np*T) by 1
    SimpiN = exp(rev_intcpt +(1-sigmahat)*(betak*Sim.k(find(Sim.year==(t))) - Sim.omega(find(Sim.year==(t))) )) ...
        *(1/sigmahat)*scale; % (np*T) by 1
    
    % compute transition for t+1 productivity
    tempk = cumsum(ones(Nf,1)*numkgrid) - numkgrid; % Nf x 1
    indexk = Sim.grid_kindex(find(Sim.year==t)) + tempk; % Nf x 1
    %                 grid_kindex  indexk
    %         firm 1      k1(3)      k1    (3)
    %         firm 2      k2(3)      8+k2  (11)
    %         firm 3      k3(8)      16+k3 (24)
    %
    %         EVN=
    %             EVN(k1,firm1)           EVN(k1,firm2)           EVN(k1,firm3) ,         .... , EVN(k1,firmNf)
    %             EVN(k2,firm1)           EVN(k2,firm2)           EVN(k2,firm3) ,         .... , EVN(k2,firmNf)
    %             EVN(k3,firm1)<<<(3)     EVN(k3,firm2) <<<(11)   EVN(k3,firm3) ,         .... , EVN(k2,firmNf)
    %                 .
    %                 .
    %                 .
    %             EVN(k8,firm1)           EVN(k8,firm2)           EVN(k8,firm3)<<<(24) ,  .... , EVN(k8,firmNf)
    
    max_xgrid = find(xgrid==max(xgrid));
    min_xgrid = find(xgrid==min(xgrid));
    
    % based on period t productivity
    temp = repmat(xgrid',Nf,1) - alpha0 ...
        - alpha1*repmat(Sim.omega(find(Sim.year==(t))),1,numprodgrid)...
        - alpha2*repmat(Sim.omega(find(Sim.year==(t))).^2,1,numprodgrid)...
        - alpha3*repmat(Sim.omega(find(Sim.year==(t))).^3,1,numprodgrid); % Nf by ngrid
    
    
    % Transition probability for nonimporter
    pxN = normpdf((temp),0,sigmaxi);          % (np*T) by ngrid
    TxNobs= pxN./repmat(sum(pxN,2),1,numprodgrid);
    
    low  = find( sum(pxN,2)==0 & sign(temp(:,1))< 0);
    high = find( sum(pxN,2)==0 & sign(temp(:,1))> 0);
    
    TxNobs(find(sum(pxN,2)==0),:) = 0;
    TxNobs(low,min_xgrid) = 1;
    TxNobs(high,max_xgrid) = 1;
    clear low high
    
    % Transition probability for importer
    pxI = normpdf((temp-alpha4),0,sigmaxi);          % Nf x 100
    TxIobs= pxI./repmat(sum(pxI,2),1,numprodgrid);
    
    low  = find( sum(pxI,2)==0 & sign(temp(:,1))< 0);
    high = find( sum(pxI,2)==0 & sign(temp(:,1))> 0);
    
    TxIobs(find(sum(pxN,2)==0),:) = 0;
    TxIobs(low,min_xgrid) = 1;
    TxIobs(high,max_xgrid) = 1;
    
    clear temp
    
    % Expected continuation value
    SimEVN = (TxNobs* VN )'; % (Nfx100 * 100x8)' = 8xNf
    SimEVN = SimEVN(indexk); % nFx1
    
    SimEVI = (TxIobs* VI)'; % nkgrid by (np*T)
    SimEVI = SimEVI(indexk); % (np*T) by 1
    
    % Marginal benefit
    Sim.dIN  = SimpiI - SimpiN + delta*(SimEVI - SimEVN); % (np*T) by 1
    Sim.sign = sign(Sim.dIN-Sim.cost);
    
    % decide import or not (import if )
    Sim.imptemp = ones(Nf,1);
    Sim.imptemp(find(Sim.sign<0)) = 0;
    
    Sim.imp(find(Sim.year==t)) = Sim.imptemp;
    disp(['Period' num2str(t,'%10.0f') 'simulation completed'])
end
    



%% Export participation

% data
ImpPar.data = ones(9,Sim.T);
for t = 1:Sim.T
y = t + 1999;
ImpPar.data(1,t) = mean(imp_dummy(find(year==y)));
ImpPar.data(2,t) = mean(imp_dummy(find(year==y & grid_kindex == 1)));
ImpPar.data(3,t) = mean(imp_dummy(find(year==y & grid_kindex == 2)));
ImpPar.data(4,t) = mean(imp_dummy(find(year==y & grid_kindex == 3)));
ImpPar.data(5,t) = mean(imp_dummy(find(year==y & grid_kindex == 4)));
ImpPar.data(6,t) = mean(imp_dummy(find(year==y & grid_kindex == 5)));
ImpPar.data(7,t) = mean(imp_dummy(find(year==y & grid_kindex == 6)));
ImpPar.data(8,t) = mean(imp_dummy(find(year==y & grid_kindex == 7)));
ImpPar.data(9,t) = mean(imp_dummy(find(year==y & grid_kindex == 8)));
end

% data
ImpPar.sim = ones(9,Sim.T);
for t = 1:Sim.T
y = t;
ImpPar.sim(1,t) = mean(Sim.imp(find(Sim.year==y)));
ImpPar.sim(2,t) = mean(Sim.imp(find(Sim.year==y & Sim.grid_kindex == 1)));
ImpPar.sim(3,t) = mean(Sim.imp(find(Sim.year==y & Sim.grid_kindex == 2)));
ImpPar.sim(4,t) = mean(Sim.imp(find(Sim.year==y & Sim.grid_kindex == 3)));
ImpPar.sim(5,t) = mean(Sim.imp(find(Sim.year==y & Sim.grid_kindex == 4)));
ImpPar.sim(6,t) = mean(Sim.imp(find(Sim.year==y & Sim.grid_kindex == 5)));
ImpPar.sim(7,t) = mean(Sim.imp(find(Sim.year==y & Sim.grid_kindex == 6)));
ImpPar.sim(8,t) = mean(Sim.imp(find(Sim.year==y & Sim.grid_kindex == 7)));
ImpPar.sim(9,t) = mean(Sim.imp(find(Sim.year==y & Sim.grid_kindex == 8)));
end

r = 150; % pixels per inch
set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [0 0 1280 1024]/r);
print(gcf,'-dpng',sprintf('-r%d',r), 'F01ExpProdAut.png');

graphymin = 2000
graphymax = 2000 + Sim.T -1;
% Compare result
h=plot(graphymin:graphymax,ImpPar.data(1,:),'k',graphymin:graphymax,ImpPar.sim(1,:),'k:' )
title('Import Participation Rate (all firms)', 'FontSize', 20)
xlabel('Year','FontSize', 16)
ylabel('Import Participation rate','FontSize', 16)
legend({'Data','Model'},'Location','southeast','fontsize',16)
set(gca,'FontSize',16);
set(h,{'LineWidth'},{4;4})
print(gcf,'-dpng',sprintf('-r%d',r), 'FigDiff00.png');

BaseName='FigDiff0';
for k=1:8
s = k+1
FileName=[BaseName,num2str(k)]
h=plot(graphymin:graphymax,ImpPar.data(s,:),'k',graphymin:graphymax,ImpPar.sim(s,:),'k:' )
title(['Import Participation Rate (',num2str(k),') '], 'FontSize', 20)
xlabel('Year','FontSize', 16)
ylabel('Import Participation rate','FontSize', 16)
legend({'Data','Model'},'Location','southeast','fontsize',16)
set(gca,'FontSize',16);
set(h,{'LineWidth'},{4;4})
print(gcf,'-dpng',sprintf('-r%d',r), FileName);
end
