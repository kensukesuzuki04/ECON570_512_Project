% Nonimporter profit
%run B05_NonimporterProfit.m

global piI numprodgrid kgrid ngrid nhsgrid betak betan xgrid rev_intcpt sigmahat 


piI = zeros(numprodgrid,length(kgrid)); % 
temp1=betak*repmat(kgrid',numprodgrid,1)  ...
        + betan*repmat(ngrid',numprodgrid,1) ...  ---> importer 1
        - repmat(xgrid,1,8);
piI=exp(rev_intcpt+(1-sigmahat)*temp1)*(1/sigmahat)*scale;

% piI = zeros(numprodgrid,length(kgrid)); % 
% temp1=betak*repmat(kgrid',numprodgrid,1)  ...
%         + betan*repmat(ngrid',numprodgrid,1) ...  ---> importer 1
%         - repmat(xgrid,1,8);
% piI=exp(rev_intcpt+(1-sigmahat)*temp1)*(1/sigmahat)*scale;

%pih=exp(hrcoast+(1+etah)*temp1)*(-1/etah)*scale;

clear temp1
%EOF

