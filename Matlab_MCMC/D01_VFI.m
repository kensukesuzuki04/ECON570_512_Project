% Value Function Iteration

%

%Choice specific continuation values
EVN = zeros(numprodgrid,numkgrid); % nonimporter
EVI = zeros(numprodgrid,numkgrid); % importer

%Choice probabilities
PNN = zeros(numprodgrid,numkgrid); % Non --> Non
PNI = zeros(numprodgrid,numkgrid); % Non --> Imp
PIN = zeros(numprodgrid,numkgrid); % Imp --> Non
PII = zeros(numprodgrid,numkgrid); % Imp --> Imp

%Initialize the values
VN = zeros(numprodgrid,numkgrid); % nonimporter
VI = zeros(numprodgrid,numkgrid); % importer

%Starting Values 
VNnew = piN/(1-delta);
VInew = (piI)/(1-delta);


% %Define the cost vectors
% gammaIF = 2;
% marginS = 30;
% gammaNS =  gammaIF + marginS;

iter=0;


while max(max(abs(VN-VNnew)))>1e-6||...
      max(max(abs(VI-VInew)))>1e-6

VN = VNnew;  
VI = VInew;  

% expected continuation value (choose Nonimporter)
EVH = TxN*VN; 
% expected continuation value (choose Importer)
EVI = TxI*VI; 

% marginal benefit
dIN   = piI - piN + delta*(EVI - EVN);
%dIN(find(dIN<0)) = 1e-6;

%CHOICE PROBABILITIES
% (1) Starting as nonimporter
PNN = F03condProb(dIN, gammaNS);
PNI = 1 - PNN;

% (2) Start as nonimporter
PIN = F03condProb(dIN,gammaIF);
PII = 1 - PIN;

% interate on UC.VH, UC.VI, UC.VD, C.VH, C.VI
% (1) value of starting as nonimporter
%condENI= gammaNS*ones(100,8) - dIN.*exp(-dIN/gammaNS)./expcdf(dIN,gammaNS);
%*gammaNS;
%condENI_ = expcdf(dIN,gammaNS)*gammaNS - dIN.*exp(-1/gammaNS*dIN)./expcdf(dIN,gammaNS);
% VNnew = PNN.*(piN + delta*EVN)+ ...
%         PNI.*(piI_mean - condENI + delta*EVI);
VNnew = piI + delta*EVI - expcdf(dIN,gammaNS)*gammaNS;
%VNnew - VNnew_
    
% (1) value of starting as importer
%condEII=expcdf(dIN,gammaIF)*gammaIF;
%condEII= gammaIF*ones(100,8) - dIN.*exp(-dIN/gammaIF)./expcdf(dIN,gammaIF);
VInew = piI + delta*EVI - expcdf(dIN,gammaIF)*gammaIF;

diff_x = max([ max(max(abs(VN-VNnew))); max(max(abs(VI-VInew)))]);
iter = iter+1;
end

%fprintf('VFI converged. Time epalsed was %d', valuefun_time)
%disp('')


%EOF
