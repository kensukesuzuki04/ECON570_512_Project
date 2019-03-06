%run B02_lnk_lnn_grid.m

global xgrid omega numprodgrid

% Descretize Productivity --- Random grid
u = haltonseq(1,numprodgrid)';

xmin = min(omega);
xmax = max(omega);
xgrid = xmin + (xmax - xmin)*u; % 100 by 1