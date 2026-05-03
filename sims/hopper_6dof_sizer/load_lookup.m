mdot_data = readmatrix("mdot_lookup.xlsx");

thrust_data = mdot_data(:,1);   % breakpoints
fu_mdot_data = mdot_data(:,2);
ox_mdot_data = mdot_data(:,3);

%load wind
load('wind_vectors2.mat');

% Load CG MOI Tables
load("cg_I_LUT.mat");