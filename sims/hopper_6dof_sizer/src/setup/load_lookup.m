mdot_data = readmatrix("mdot_lookup.xlsx");

thrust_data = mdot_data(:,1);   % breakpoints
fu_mdot_data = mdot_data(:,2);
ox_mdot_data = mdot_data(:,3);

%load wind
load('wind_vectors2.mat');

% Load CG MOI Tables
load("cg_I_LUT.mat");

% --- Load Engine Performance Data ---
mdot_data = readmatrix("mdot_lookup.xlsx");
thrust_data = mdot_data(:,1);   % breakpoints
fu_mdot_data = mdot_data(:,2);
ox_mdot_data = mdot_data(:,3);

assignin('base', 'thrust_data',  thrust_data);
assignin('base', 'fu_mdot_data', fu_mdot_data);
assignin('base', 'ox_mdot_data', ox_mdot_data);


% --- Load Wind Data Dynamically ---
wind_data = load('wind_vectors2.mat');
wind_fields = fieldnames(wind_data);
for i = 1:length(wind_fields)
    assignin('base', wind_fields{i}, wind_data.(wind_fields{i}));
end


% --- Load CG & MOI Tables Dynamically ---
cg_data = load('cg_I_LUT.mat');
cg_fields = fieldnames(cg_data);
for i = 1:length(cg_fields)
    assignin('base', cg_fields{i}, cg_data.(cg_fields{i}));
end