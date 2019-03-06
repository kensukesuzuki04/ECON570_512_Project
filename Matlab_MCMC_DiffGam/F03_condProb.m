function funph=F03condProb(grid_kindex_, MB,gamma)

global numkgrid

% draw different fixed/sunk cost depending on capital

funph = zeros(size(MB));

for cap = 1:numkgrid % for each 
    gammacap = gamma(cap); % take relevant gamma
    tempindex = find(grid_kindex_ == cap);
    funph(tempindex) = (1-expcdf(MB(tempindex),gammacap));
end

end
%End of File
