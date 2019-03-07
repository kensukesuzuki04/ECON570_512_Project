% obtain standard error

clear all

load R_MLEpara

% Use HS
usehs = 0;
usegridn = 0;

% ISIC code
isic = 29;

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

F04_loglf(fmin_parameter)

jacobian(gradient(F04_loglf))