%run B01_FirstStage.m

% Declar global
global mean_lnk1 kgrid ngrid grid_lnk1 grid_lnn1  grid_kindex grid_kindex_ grid_lnn1_ ...
    grid_lnnhs1 grid_lnnhs1_ grid_lnn1_cand grid_lnnhs1_cand nhsgrid grid_lnk1_

%% capital grid

% Pick quantile (12.5%)
quantilelnk1 = quantile(lnk1,[0 0.125 0.25 0.375 0.5 0.625 0.75 0.875 1]);

size=zeros(8,1);
size(1) = length(find(lnk1>=quantilelnk1(1) & lnk1<quantilelnk1(2)));
size(2) = length(find(lnk1>=quantilelnk1(2) & lnk1<quantilelnk1(3) ));
size(3) = length(find(lnk1>=quantilelnk1(3) & lnk1<quantilelnk1(4) ));
size(4) = length(find(lnk1>=quantilelnk1(4) & lnk1<quantilelnk1(5) ));
size(5) = length(find(lnk1>=quantilelnk1(5) & lnk1<quantilelnk1(6) ));
size(6) = length(find(lnk1>=quantilelnk1(6) & lnk1<quantilelnk1(7) ));
size(7) = length(find(lnk1>=quantilelnk1(7) & lnk1<quantilelnk1(8) ));
size(8) = length(find(lnk1>=quantilelnk1(8)));

disp ('--- Number of observations in each bin ---')
size

% Make grid: get mean of each bin
grid=zeros(8,1);
grid(1) = mean(lnk1(find(lnk1>=quantilelnk1(1) & lnk1<quantilelnk1(2) )));
grid(2) = mean(lnk1(find(lnk1>=quantilelnk1(2) & lnk1<quantilelnk1(3) )));
grid(3) = mean(lnk1(find(lnk1>=quantilelnk1(3) & lnk1<quantilelnk1(4) )));
grid(4) = mean(lnk1(find(lnk1>=quantilelnk1(4) & lnk1<quantilelnk1(5) )));
grid(5) = mean(lnk1(find(lnk1>=quantilelnk1(5) & lnk1<quantilelnk1(6) )));
grid(6) = mean(lnk1(find(lnk1>=quantilelnk1(6) & lnk1<quantilelnk1(7) )));
grid(7) = mean(lnk1(find(lnk1>=quantilelnk1(7) & lnk1<quantilelnk1(8) )));
grid(8) = mean(lnk1(find(lnk1>=quantilelnk1(8) & lnk1<=quantilelnk1(9) )));

% grid point
kgrid = grid;

% Compute mean of log k for each firm
firm_id_list = unique(id);
mean_lnk1 = zeros(N,1);
for i = 1:Nf
    j = firm_id_list(i);
    mean_lnk1(find(id==j)) = mean(lnk1(find(id==j)));
end

grid_lnk1=zeros(N,1);
grid_lnk1(find(mean_lnk1>=quantilelnk1(1) & mean_lnk1<quantilelnk1(2) )) = grid(1);
grid_lnk1(find(mean_lnk1>=quantilelnk1(2) & mean_lnk1<quantilelnk1(3) )) = grid(2);
grid_lnk1(find(mean_lnk1>=quantilelnk1(3) & mean_lnk1<quantilelnk1(4) )) = grid(3);
grid_lnk1(find(mean_lnk1>=quantilelnk1(4) & mean_lnk1<quantilelnk1(5) )) = grid(4);
grid_lnk1(find(mean_lnk1>=quantilelnk1(5) & mean_lnk1<quantilelnk1(6) )) = grid(5);
grid_lnk1(find(mean_lnk1>=quantilelnk1(6) & mean_lnk1<quantilelnk1(7) )) = grid(6);
grid_lnk1(find(mean_lnk1>=quantilelnk1(7) & mean_lnk1<quantilelnk1(8) )) = grid(7);
grid_lnk1(find(mean_lnk1>=quantilelnk1(8) & mean_lnk1<=quantilelnk1(9) )) = grid(8);

% it may vary across years
grid_lnk1_=zeros(N,1);
grid_lnk1_(find(lnk1>=quantilelnk1(1) & lnk1<quantilelnk1(2) )) = grid(1);
grid_lnk1_(find(lnk1>=quantilelnk1(2) & lnk1<quantilelnk1(3) )) = grid(2);
grid_lnk1_(find(lnk1>=quantilelnk1(3) & lnk1<quantilelnk1(4) )) = grid(3);
grid_lnk1_(find(lnk1>=quantilelnk1(4) & lnk1<quantilelnk1(5) )) = grid(4);
grid_lnk1_(find(lnk1>=quantilelnk1(5) & lnk1<quantilelnk1(6) )) = grid(5);
grid_lnk1_(find(lnk1>=quantilelnk1(6) & lnk1<quantilelnk1(7) )) = grid(6);
grid_lnk1_(find(lnk1>=quantilelnk1(7) & lnk1<quantilelnk1(8) )) = grid(7);
grid_lnk1_(find(lnk1>=quantilelnk1(8) & lnk1<=quantilelnk1(9) )) = grid(8);

% constant grid
grid_kindex=zeros(N,1);
grid_kindex(find(mean_lnk1>=quantilelnk1(1) & mean_lnk1<quantilelnk1(2) )) = 1;
grid_kindex(find(mean_lnk1>=quantilelnk1(2) & mean_lnk1<quantilelnk1(3) )) = 2;
grid_kindex(find(mean_lnk1>=quantilelnk1(3) & mean_lnk1<quantilelnk1(4) )) = 3;
grid_kindex(find(mean_lnk1>=quantilelnk1(4) & mean_lnk1<quantilelnk1(5) )) = 4;
grid_kindex(find(mean_lnk1>=quantilelnk1(5) & mean_lnk1<quantilelnk1(6) )) = 5;
grid_kindex(find(mean_lnk1>=quantilelnk1(6) & mean_lnk1<quantilelnk1(7) )) = 6;
grid_kindex(find(mean_lnk1>=quantilelnk1(7) & mean_lnk1<quantilelnk1(8) )) = 7;
grid_kindex(find(mean_lnk1>=quantilelnk1(8) & mean_lnk1<=quantilelnk1(9))) = 8;

% index -- may vary over time
grid_kindex_=zeros(N,1);
for i = 1:8
grid_kindex_(find(grid_lnk1_ == grid(i) )) = i;
end

clear grid size quantile

%% n grid

ngrid = zeros(8,1);
for i = 1:8
       
    ngrid(i) = median(lnn1(find(grid_lnk1_==kgrid(i) & imp_dummy ==1))) % importer with capital grid 1
    
end

grid_lnn1_ = zeros(N,1);
grid_lnn1 = zeros(N,1);
grid_lnn1_cand = zeros(N,1);

for i = 1:8
    tempk = kgrid(i);
    grid_lnn1_(find(grid_lnk1_ == tempk & imp_dummy ==1 )) =  ngrid(i);
    
    grid_lnn1(find(grid_lnk1 == tempk  & imp_dummy ==1)) =  ngrid(i);
    grid_lnn1_cand(find(grid_lnk1 == tempk)) =  ngrid(i);
end

%% use hs
nhsgrid = zeros(8,1);
for i = 1:8
    nhsgrid(i) = mean(lnnhs1(find(grid_lnk1_==kgrid(i) & imp_dummy ==1))); % importer with capital grid 1
end

grid_lnnhs1_ = zeros(N,1);
grid_lnnhs1 = zeros(N,1);
grid_lnnhs1_cand = zeros(N,1);

for i = 1:8
    tempk = kgrid(i);
    grid_lnnhs1_(find(grid_lnk1_ == tempk & imp_dummy ==1 )) =  nhsgrid(i);
    grid_lnnhs1(find(grid_lnk1 == tempk  & imp_dummy ==1)) =  nhsgrid(i);
    grid_lnnhs1_cand(find(grid_lnk1 == tempk)) =   nhsgrid(i);
end
    

clear grid size quantile tempk


