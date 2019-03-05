% Function returns the residual

function fval = F01_NLSevol_IV (para, lagphi, laglnk, laglnn, lagimp_dummy, lnk, lnn, phi, sigma, NwoF)

% note that NwoF is number of obervation excluding the first year 
% ("N"umber of obervation "w/o" "F"irst year)
% Since we use lagged variables, samples are restricted to non-first year

alpha0 = para(1,1);   % extract alpha0 - alpha 4 
alpha1 = para(2,1);
alpha2 = para(3,1);
alpha3 = para(4,1);
alpha4 = para(5,1);
betak = para(length(para)-1,1);      % etak
betan = para(length(para),1);        % etan

betaks = (1-sigma)*betak;
betans = (1-sigma)*betan;
alpha0s = -(1-sigma)*alpha0;
alpha2s = - alpha2/(1-sigma);
alpha3s = alpha3/((1-sigma)^2);
alpha4s = -(1-sigma)*alpha4;
alphas = [alpha0s; alpha1; alpha2s; alpha3s; alpha4s];

lagphikn = lagphi - betaks * laglnk - betans * laglnn ;
lagphikn2 = lagphikn.^2;
lagphikn3 = lagphikn.^3;


RHS_X = [ones(NwoF,1), lagphikn, lagphikn2, lagphikn3, lagimp_dummy];


RHS = betaks * lnk  + betans* lnn + RHS_X*alphas;
LHS = phi ;

fval = RHS-LHS;


end