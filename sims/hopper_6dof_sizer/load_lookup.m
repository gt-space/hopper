data = readmatrix("mdot_lookup.xlsx");

bp = data(:,1);   % breakpoints
tbl1 = data(:,2);
tbl2 = data(:,3);

%load wind
load('wind_vectors2.mat');

% Load CG MOI Tables
load("cg_I_LUT.mat");