
clearvars;

main();

load_lookup();

cg_moi_test();

LinerizationMaster();

% Baffle damping LUT: rerun generateBaffleLUT after changing tanks/baffles
if ~isfile('baffle_lut.csv')
    generateBaffleLUT();
end
baffle_lut = readtable('baffle_lut.csv');
for lut_var = baffle_lut.Properties.VariableNames
    assignin('base', lut_var{1}, baffle_lut.(lut_var{1})');
end

if abs(max(ox_mass_profile) - IN.propulsion.oxidizer_mass) > 1e-6 || ...
   abs(max(fu_mass_profile) - IN.propulsion.fuel_mass) > 1e-6
    warning('baffle_lut.csv was generated for different propellant masses; rerun generateBaffleLUT.');
end

%prop_system_testing();

    