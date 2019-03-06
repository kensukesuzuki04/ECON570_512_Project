function out = F04_loglf(par)
global Nf S T scale  notfirst lastyear numprodgrid ...
    xgrid lnk1 grid_lnk1 grid_lnn1 grid_lnnhs1 grid_lnn1_cand ...
    grid_kindex grid_kindex_   omega ...
    piN piI TxN TxI ...
    imp_dummy lagimp_dummy ...
    numkgrid ...
    sigmaxi sigmahat delta ...
    betak betan alpha0 alpha1 alpha2 alpha3 alpha4 ...
    TxIobs TxNobs ...
    rev_intcpt 


%v1parameters;

%Dynamic parameters --- par = [gammaNS; gammaIF]
%Costs faced by Non-Exporter

gammaIF = par(1:8);
gammaNS = par(9:16);
a0 = par(17);
a1 = par(18);
a2 = par(19);
%marginS=par(2);

% gammaIF = 2;
% marginS = 30;
%gammaNS =  gammaIF + marginS;

% %Export revenue intercept

%D Valuefunction Iteration
D01_VFI;

tempk = cumsum(ones(Nf*T,1)*numkgrid)-numkgrid; % (np*T) by 1
indexk = grid_kindex+tempk; % (np*T) by 1

% forward lag grid capital stock
lagindexk_ = zeros(T*Nf,1);
lagindexk_(1:T*Nf-1) = grid_kindex_(2:T*Nf);
lagindexk_ = lagindexk_ .* (1-lastyear);
% 2007 expected capital stock is same as 2006
lagindexk_(find(lastyear == 1)) = grid_kindex_(find(lastyear==1));

tempk = cumsum(ones(Nf*T,1)*numkgrid)-numkgrid; % (np*T) by 1
indexk_ = lagindexk_+tempk; % (np*T) by 1

    max_xgrid = find(xgrid==max(xgrid));
    min_xgrid = find(xgrid==min(xgrid));

lf= ones(Nf,S);



for s = 1:1:S
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %1.DATA
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    x  = omega; % (np*T) by 1 ... this is data
    k  = grid_lnk1;    % (np*T) by 1 ... grid point
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %2. OBSERVATION-TIME SPECIFIC p(x_n|x_it,m_it)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % p(x_n|x_it,d_it)
    % given productivity estimates (data), get transition probability
    
    temp = repmat(xgrid',Nf*T,1) - alpha0 ...
        - alpha1*repmat(x,1,numprodgrid)...
        - alpha2*repmat(x.^2,1,numprodgrid)...
        - alpha3*repmat(x.^3,1,numprodgrid); % (np*T) by ngrid
    
    % Nonimporter
    pxN = normpdf(temp,0,sigmaxi);          % (np*T) by ngrid
    tempden = repmat(sum(pxN,2),1,numprodgrid);
    TxNobs= pxN./repmat(sum(pxN,2),1,numprodgrid);
    
    low  = find( sum(pxN,2)==0 & sign(temp(:,1))< 0);
    high = find( sum(pxN,2)==0 & sign(temp(:,1))> 0);

    TxNobs(find(sum(pxN,2)==0),:) = 0;
    TxNobs(low,min_xgrid) = 1;
    TxNobs(high,max_xgrid) = 1;
    clear low high

    % Importer
    temp1 = temp-alpha4;
    pxI = normpdf((temp1),0,sigmaxi);          % (np*T) by ngrid
    TxIobs= pxI./repmat(sum(pxI,2),1,numprodgrid);

    low  = find( sum(pxI,2)==0 & sign(temp1(:,1))< 0);
    high = find( sum(pxI,2)==0 & sign(temp1(:,1))> 0);

    TxIobs(find(sum(pxN,2)==0),:) = 0;
    TxIobs(low,min_xgrid) = 1;
    TxIobs(high,max_xgrid) = 1;

    clear low high temp1
    
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %3. CALCULATE EV_it+1 pixi_t pixd_t dIHit+1 dDHit+1 dIDit+1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% EXPECTED CONTINUATION VALUE
    % Use value function solved above
    OBS.EVN = (TxNobs* VN )'; % nkgrid by (np*T)
    % use actual tomorrow's capital stock
    OBS.EVN = OBS.EVN(indexk); % (np*T) by 1
    
    OBS.EVI = (TxIobs* VI)'; % nkgrid by (np*T)
    OBS.EVI = OBS.EVI(indexk); % (np*T) by 1
    
    % PROFIT
    % Importer use mean_lnn1 for cost reduction effect
    piIobs = exp(rev_intcpt +(1-sigmahat)*(betak*grid_lnk1 + betan * grid_lnn1_cand - x ))*(1/sigmahat)*scale; % (np*T) by 1
    piNobs = exp(rev_intcpt +(1-sigmahat)*(betak*grid_lnk1 - x ))*(1/sigmahat)*scale; % (np*T) by 1
    
    % Marginal benefit
    OBS.dIN  = piIobs - piNobs + delta*(OBS.EVI - OBS.EVN); % (np*T) by 1
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%grid%%%%%
    %4 CHOICE PROBABILITIES
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % lnP(dit|dit-1)
    %i. start as Nonimporter
    OBS.PNN= F03_condProb(grid_kindex_, OBS.dIN,gammaNS);
    OBS.PNI = 1 - OBS.PNN;
    
    %ii. start as Importer
    OBS.PIN=F03_condProb(grid_kindex_,OBS.dIN, gammaIF);
    OBS.PII = 1 - OBS.PIN;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %5 First calculate for all NON-INITIAL years
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Unconditional Choice Probabilities
    % to become Nonimporter
    ccpN   =(lagimp_dummy==0).* OBS.PNN + (lagimp_dummy==1).*OBS.PIN;
    
    % to become Importer
    ccpI   =(lagimp_dummy==0).* OBS.PNI + (lagimp_dummy==1).*OBS.PII;
    
    lmnonini =notfirst.* ((imp_dummy==0).*ccpN + (imp_dummy==1).*ccpI );
    lmnonininew = reshape(lmnonini,T,Nf);
    lmnonininew = lmnonininew(2:T,:);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %6 Second calculate for INITIAL years
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %5.1 (di0==H).*lnP(di0=H|.)+(di0==I).*lnP(di0==I|.)+(di0==D).*lnP(di0==D|.);
    %proxy P(di0==H|.) using a simple logit
    
    OBS.summ = 1+exp(a0 + a1*x + a2*grid_lnk1);
    OBS.P0N = 1./OBS.summ;
    OBS.P0I = exp(a0 + a1*x + a2*grid_lnk1)./OBS.summ;
    
    lmini= (1-notfirst).* ((imp_dummy==0).*OBS.P0N + (imp_dummy==1).*OBS.P0I);
    lmininew = reshape(lmini,T,Nf);
    lmininew = lmininew(1,:); 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %7 calculate p(mi0^T|xi0^T,ki0^T)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    lmnew = [lmininew;lmnonininew];
    
%     lmnew = reshape(lm,T,Nf);
%     lmnew = lmnew(2:T,:);
    li = (prod(lmnew))';
    

    lf(:,s) = max(li,1e-30);
end
ml = lf;
% sgm = mean(gm,2);
% sgz = mean(gz,2);
% outm = -sum(log(sgm));
% outz = -sum(log(sgz));
% [outm outz];

out = -sum(log(ml));









