%run B03_RandomGrid.m

%Transitions of productivity x on gridx
global TxI TxN alpha0 alpha1 alpha2 alpha3 alpha4 numprodgrid sigmaxi % importer and nonimporter

TxI = zeros(numprodgrid,numprodgrid); % 100 x 100 for importer
TxN = zeros(numprodgrid,numprodgrid); % 100 x 100 for nonimporter
 
temp1=repmat(xgrid,1,numprodgrid); 
temp2 =repmat(xgrid',numprodgrid,1) ...
        -alpha0 ...
        -alpha1*temp1...
        -alpha2*temp1.^2 ...
        -alpha3*temp1.^3; % xi for nonimporters
 
% Nonimporter
TxN= normpdf(temp2,0,sigmaxi); % 100 x 100 f(xi_ij|omega_i, Nonimporter)
TxN= TxN./repmat(sum(TxN,2),1,numprodgrid);

% Importer
TxI= normpdf((temp2-alpha4),0,sigmaxi);
TxI= TxI./repmat(sum(TxI,2),1,numprodgrid);

% Tx0 = zeros(numprodgrid,numprodgrid);
% Tx1 = zeros(numprodgrid,numprodgrid);
% for i=1:1:numprodgrid
%     temp=xgrid-Para.alpha0-Para.alpha1*xgrid(i)-Para.alpha2*xgrid(i)^2-Para.alpha3*xgrid(i)^3;
%     Tx0(i,:)=normpdf((temp)/Para.sigmaxi)'/Para.sigmaxi;
%     Tx0(i,:)=Tx0(i,:)/sum(Tx0(i,:));
%     Tx1(i,:)=normpdf((temp-Para.alpha4)/Para.sigmaxi)'/Para.sigmaxi;
%     Tx1(i,:)=Tx1(i,:)/sum(Tx1(i,:));
% end
% 
% diffI =  Tx0 - TxN;

clear temp1 temp2
