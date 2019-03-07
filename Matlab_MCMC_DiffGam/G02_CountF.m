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
Sim.T = 30;

useMCMC = 1;
usehs = 0;
usegridn = 1;
isic = 29;

% read parameter either from MCMC or MLE
if useMCMC == 1
    Sim.gammaIF = ps_mean(1:8);
    Sim.gammaNS = ps_mean(9:16);
else
    Sim.gammaIF = fmiunc.par(1);
    Sim.gammaNS = fminunc.par(2);
end

counterfactual = 1;
% 1 Uniform FS 10%
% 2 Uniform F 10%
% 3 Uniform S 10%
% 4 Small FS 10%
% 5 Small F 10%
% 6 Small S 10%

% subsidy in percentage
sub = 0.20;
subsmall = 0.20;


%% Go over First Stage Estimation

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


%% setup for simulation

Sim.year = repmat([1:1:Sim.T]',Nf,1); 
Sim.k = zeros(Nf*Sim.T, 1)
Sim.kind = zeros(Nf*Sim.T, 1)
Sim.n_cand = zeros(Nf*Sim.T, 1)
Sim.omega = zeros(Nf*Sim.T, 1) % factual 
Sim.omegaC = zeros(Nf*Sim.T, 6) % counterfactual
Sim.imp = 999*ones(Nf*Sim.T, 1) % factual
Sim.impC = 999*ones(Nf*Sim.T, 6) % counterfactual
Sim.grid_kindex = zeros(Nf*Sim.T, 1)
Sim.FYgrid_kindex = zeros(Nf,1)

Sim.subC = zeros(Nf*Sim.T,6)% counterfactual

% assign capital (grid) and candidate n (grid) from data
for i = 1:Nf
    Sim.k((i-1)*Sim.T+1:i*Sim.T) = grid_lnk1((i-1)*T+1,1);
    Sim.grid_kindex((i-1)*Sim.T+1:i*Sim.T) = grid_kindex((i-1)*T+1,1);
    Sim.FYgrid_kindex = Sim.grid_kindex(find(Sim.year==1));
    Sim.n_cand((i-1)*Sim.T+1:i*Sim.T) = grid_lnn1_cand((i-1)*T+1,1); % grid_lnn1_cand assumes constant capital stock
end

% assing first year productivity
Sim.omega(find(Sim.year==1)) = omega(find(year==2000));
Sim.omegaC(find(Sim.year==1),:) = repmat(omega(find(year==2000)),1,6);

% use first year import status
Sim.imp(find(Sim.year==1)) = imp_dummy(find(year==2000));
Sim.impC(find(Sim.year==1),:) = repmat(imp_dummy(find(year==2000)),1,6);


%% Simulation with estimates

gammaIF = Sim.gammaIF;
gammaNS = Sim.gammaNS;
D01_VFI


% compute next period productivity (based on today's productivity and import status)
for t = 2:Sim.T
    
    seed= 19900904*t;
    rng(seed);
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
    gammaIFCcap = Sim.gammaIF(cap);
    gammaNSCcap = Sim.gammaNS(cap);
    Sim.cost(find(Sim.FYgrid_kindex==cap)) =  Sim.imp(find(Sim.year==(t-1) & Sim.grid_kindex==cap)) .* exprnd(gammaIFCcap,[length(Sim.imp(find(Sim.year==(t-1) & Sim.grid_kindex==cap))),1]) ... % case of importer
        + (1-Sim.imp(find(Sim.year==(t-1) & Sim.grid_kindex==cap))) .* exprnd(gammaNSCcap,[length(Sim.imp(find(Sim.year==(t-1) & Sim.grid_kindex==cap))),1]); % case of  non importer
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

%% Simulation with counterfactual

for counterfactual = 1:6
    
    % counter factual parameters
    Sim.CgammaIF = Sim.gammaIF;
    Sim.CgammaNS = Sim.gammaNS;
    
    % 1 Uniform FS 10%
    % 2 Uniform F 10%
    % 3 Uniform S 10%
    % 4 Small FS 10%
    % 5 Small F 10%
    % 6 Small S 10%
    
    if counterfactual == 1
        Sim.CgammaIF = (1-sub) .* Sim.CgammaIF;
        Sim.CgammaNS = (1-sub) .* Sim.CgammaNS;
    elseif counterfactual == 2
        Sim.CgammaIF = (1-sub) .*Sim.CgammaIF;
    elseif counterfactual == 3
        Sim.CgammaNS = (1-sub) .*Sim.CgammaNS;
    elseif counterfactual == 4
        Sim.CgammaIF(1) = (1-subsmall) .* Sim.CgammaIF(1);
        Sim.CgammaNS(1) = (1-subsmall) .* Sim.CgammaNS(1);
    elseif counterfactual == 5
        Sim.CgammaIF(1) = (1-subsmall) .* Sim.CgammaIF(1);
    else
        Sim.CgammaNS(1) = (1-subsmall) .* Sim.CgammaNS(1);
    end
    
    counterfactual;
    Sim.gammaIF;
    gammaIF = Sim.CgammaIF;
    Sim.gammaNS;
    gammaNS = Sim.CgammaNS;
    
    D01_VFI
    
    
    % compute next period productivity (based on today's productivity and import status)
    for t = 2:Sim.T
        
    seed= 19900904*t;
    rng(seed);
        
        % period t productiviy
        shock = normrnd(0,sigmaxi,[Nf,1]);
        Sim.omegaC(find(Sim.year==t),counterfactual) = alpha0 + ...
            + alpha1*Sim.omegaC(find(Sim.year==(t-1)),counterfactual) ...
            + alpha2*Sim.omegaC(find(Sim.year==(t-1)),counterfactual).^2....
            + alpha3*Sim.omegaC(find(Sim.year==(t-1)),counterfactual).^3....
            + alpha4*Sim.impC(find(Sim.year==(t-1)),counterfactual) ...
            + shock;
        
        Sim.cost = zeros(Nf, 1);
        Sim.cost_if = zeros(Nf, 1);
        
        for cap = 1:numkgrid
            % draw relevant cost (depending on t-1 import status)
            gammaIFCcap = Sim.CgammaIF(cap);
            gammaNSCcap = Sim.CgammaNS(cap);
            Sim.cost(find(Sim.FYgrid_kindex==cap)) =  Sim.impC(find(Sim.year==(t-1) & Sim.grid_kindex==cap)) .* exprnd(gammaIFCcap,[length(Sim.impC(find(Sim.year==(t-1) & Sim.grid_kindex==cap))),1]) ... % case of importer
                + (1-Sim.impC(find(Sim.year==(t-1) & Sim.grid_kindex==cap),counterfactual)) .* exprnd(gammaNSCcap,[length(Sim.impC(find(Sim.year==(t-1) & Sim.grid_kindex==cap))),1]); % case of  non importer
            
            % draw relevant cost (depending on t-1 import status)
            gammaIFcap = Sim.gammaIF(cap);
            gammaNScap = Sim.gammaNS(cap);
            Sim.cost_if(find(Sim.FYgrid_kindex==cap)) =  Sim.impC(find(Sim.year==(t-1) & Sim.grid_kindex==cap)) .* exprnd(gammaIFcap,[length(Sim.impC(find(Sim.year==(t-1) & Sim.grid_kindex==cap))),1]) ... % case of importer
                + (1-Sim.impC(find(Sim.year==(t-1) & Sim.grid_kindex==cap),counterfactual)) .* exprnd(gammaNScap,[length(Sim.impC(find(Sim.year==(t-1) & Sim.grid_kindex==cap))),1]); % case of  non importer
        end
        
         Sim.subC(find(Sim.year == t),counterfactual) = Sim.cost_if - Sim.cost;
        
        % Profit at t
        SimpiI = exp(rev_intcpt +(1-sigmahat)*(betak*Sim.k(find(Sim.year==(t))) + betan * Sim.n_cand(find(Sim.year==(t))) - Sim.omegaC(find(Sim.year==(t)),counterfactual) )) ...
            *(1/sigmahat)*scale; % (np*T) by 1
        SimpiN = exp(rev_intcpt +(1-sigmahat)*(betak*Sim.k(find(Sim.year==(t))) - Sim.omegaC(find(Sim.year==(t)),counterfactual) )) ...
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
            - alpha1*repmat(Sim.omegaC(find(Sim.year==(t)),counterfactual),1,numprodgrid)...
            - alpha2*repmat(Sim.omegaC(find(Sim.year==(t)),counterfactual).^2,1,numprodgrid)...
            - alpha3*repmat(Sim.omegaC(find(Sim.year==(t)),counterfactual).^3,1,numprodgrid); % Nf by ngrid
        
        
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
        
        Sim.impC(find(Sim.year==t),counterfactual) = Sim.imptemp;
        Sim.subC(find(Sim.year == t),counterfactual) = Sim.subC(find(Sim.year == t),counterfactual) .*  Sim.impC(find(Sim.year==t),counterfactual);
        disp(['Period' num2str(t,'%10.0f') 'simulation completed'])
    end

end

%% Compute import participation rate
    
% data
ImpPar.data = ones(1,Sim.T);
for t = 1:Sim.T
y = t + 1999;
ImpPar.data(1,t) = mean(imp_dummy(find(year==y)));
end

% Simulation with estimates
ImpPar.sim = ones(1,Sim.T);
for t = 1:Sim.T
y = t;
ImpPar.sim(1,t) = mean(Sim.imp(find(Sim.year==y)));
end

% Simulation with counterfactual
ImpPar.simC = ones(6,Sim.T);
ImpPar.subC = ones(6,Sim.T);
for t = 1:Sim.T
y = t;
for c = 1:6
    ImpPar.simC(c,t) = mean(Sim.impC(find(Sim.year==y),c));
    ImpPar.subC(c,t) = sum(Sim.subC(find(Sim.year==y),c));
end
end

ImpPar.TotalsubC = sum(ImpPar.subC,2)/1000;

%% Generate figures

r = 150; % pixels per inch
set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [0 0 1280 1024]/r);
print(gcf,'-dpng',sprintf('-r%d',r), 'F01ExpProdAut.png');

graphymin = 2000
graphymax = 2000 + Sim.T -1;


% 1 Uniform FS 10%
% 2 Uniform F 10%
% 3 Uniform S 10%
% 4 Small FS 10%
% 5 Small F 10%
% 6 Small S 10%
clf
clear str
% compare 1 and 4
h=plot(graphymin:graphymax,ImpPar.sim(1,:),'k',graphymin:graphymax,ImpPar.simC(1,:),'k:',graphymin:graphymax,ImpPar.simC(2,:),'k--' ,graphymin:graphymax,ImpPar.simC(3,:),'k*-' )
title('Import Participation Rate: Subsidy for Fixed and Sunk Cost', 'FontSize', 20)
xlabel('Year','FontSize', 16)
ylabel('Import Participation rate','FontSize', 16)
legend({['Model'],['(1) ' num2str(sub*100) '% reduc in F&S'],['(2) ' num2str(sub*100) '% reduc in F'], ['(3) ' num2str(sub*100) '% reduc in S']},'Location','southeast','fontsize',16)
set(gca,'FontSize',16);
set(h,{'LineWidth'},{3;3;3;3})
print(gcf,'-dpng',sprintf('-r%d',r), 'FigCount00.png');

clf
clear str
% compare 1 and 4
h=plot(graphymin:graphymax,ImpPar.sim(1,:),'k',graphymin:graphymax,ImpPar.simC(1,:),'k:',graphymin:graphymax,ImpPar.simC(4,:),'k--' )
title('Import Participation Rate: Subsidy for Both F and S', 'FontSize', 20)
xlabel('Year','FontSize', 16)
ylabel('Import Participation rate','FontSize', 16)
legend({['Model'],['(1) ' num2str(sub*100) '% reduc for all firms'], ['(2) ' num2str(subsmall*100) '% reduc for Size 1 firm']},'Location','southeast','fontsize',16)
set(gca,'FontSize',16);
set(h,{'LineWidth'},{3;3;3})
ylim([0 0.5])
dim = [0.2 0.5 1 0.3];
str = {'(1) Gov Exp:' num2str(ImpPar.TotalsubC(1)) , '(2) Gov Exp' num2str(ImpPar.TotalsubC(4))};
annotation('textbox',dim,'String',str,'FitBoxToText','on');
print(gcf,'-dpng',sprintf('-r%d',r), 'FigCount14.png');

clf
clear str
% compare 2 and 5
h=plot(graphymin:graphymax,ImpPar.sim(1,:),'k',graphymin:graphymax,ImpPar.simC(2,:),'k:',graphymin:graphymax,ImpPar.simC(5,:),'k--' )
title('Import Participation Rate: Subsidy for  F', 'FontSize', 20)
xlabel('Year','FontSize', 16)
ylabel('Import Participation rate','FontSize', 16)
legend({['Model'],['(1) ' num2str(sub*100) '% reduc for all firms'], ['(2) ' num2str(subsmall*100) '% reduc for Size 1 firm']},'Location','southeast','fontsize',16)
set(gca,'FontSize',16);
set(h,{'LineWidth'},{3;3;3})
ylim([0 0.5])
dim = [0.2 0.5 1 0.3];
str = {'(1) Gov Exp:' num2str(ImpPar.TotalsubC(2)) , '(2) Gov Exp' num2str(ImpPar.TotalsubC(5))};
annotation('textbox',dim,'String',str,'FitBoxToText','on');
print(gcf,'-dpng',sprintf('-r%d',r), 'FigCount25.png');

clf
clear str
% compare 3 and 6
h=plot(graphymin:graphymax,ImpPar.sim(1,:),'k',graphymin:graphymax,ImpPar.simC(3,:),'k:',graphymin:graphymax,ImpPar.simC(6,:),'k--' )
title('Import Participation Rate: Subsidy for ', 'FontSize', 20)
xlabel('Year','FontSize', 16)
ylabel('Import Participation rate','FontSize', 16)
legend({['Model'],['(1) ' num2str(sub*100) '% reduc for all firms'], ['(2) ' num2str(subsmall*100) '% reduc for Size 1 firm']},'Location','southeast','fontsize',16)
set(gca,'FontSize',16);
set(h,{'LineWidth'},{3;3;3})
ylim([0 0.5])
dim = [0.2 0.5 1 0.3];
str = {'(1) Gov Exp:' num2str(ImpPar.TotalsubC(3)) , '(2) Gov Exp' num2str(ImpPar.TotalsubC(6))};
annotation('textbox',dim,'String',str,'FitBoxToText','on');
print(gcf,'-dpng',sprintf('-r%d',r), 'FigCount36.png');
