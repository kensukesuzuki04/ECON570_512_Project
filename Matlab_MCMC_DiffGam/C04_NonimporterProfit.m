% Nonimporter profit
%run B04_ProductivityTransition.m

global piN numprodgrid kgrid betak betan xgrid rev_intcpt sigmahat

piN = zeros(numprodgrid,length(kgrid)); % 100 x 8 matrix

% temp1=betak*repmat(kgrid',numprodgrid,1)  ...
%     - betan*ngrid(1)*ones(numprodgrid, 8) ... Non importer --> log(1)=0
%     - repmat(xgrid,1,8);
temp1=betak*repmat(kgrid',numprodgrid,1)  ...
    + betan*0*ones(numprodgrid, 8) ... Non importer --> 0
    - repmat(xgrid,1,8);

piN=exp(rev_intcpt+(1-sigmahat)*temp1)*(1/sigmahat)*scale;
%pih=exp(hrcoast+(1+etah)*temp1)*(-1/etah)*scale;

clear temp1
%EOF

