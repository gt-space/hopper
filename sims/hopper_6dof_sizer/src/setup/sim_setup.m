
clearvars;

main();

load_lookup();

cg_moi_test();

LinerizationMaster();

% Baffle damping LUT: rerun generateBaffleLUT after changing tanks/baffles
if ~isfile(hopperRoot('data', 'baffle_lut.csv'))
    generateBaffleLUT();
end
baffle_lut = readtable(hopperRoot('data', 'baffle_lut.csv'));
for lut_var = baffle_lut.Properties.VariableNames
    assignin('base', lut_var{1}, baffle_lut.(lut_var{1})');
end

burns = ~isfield(IN.propulsion, 'burns_propellant') || IN.propulsion.burns_propellant;   % false: the CTV (slosh off)
if burns && (abs(max(ox_mass_profile) - IN.propulsion.oxidizer_mass) > 1e-6 || ...
             abs(max(fu_mass_profile) - IN.propulsion.fuel_mass) > 1e-6)
    warning('baffle_lut.csv was generated for different propellant masses; rerun generateBaffleLUT.');
end
% ... and for this vehicle's tanks and baffles (the config file)
if burns
    B = load(hopperRoot('data', 'baffle_lut.mat'), 'tank_r', 'ox_tank_h', 'fu_tank_h', 'ox_baffle', 'fu_baffle');
    VT = IN.vehicle.tanks;
    built = [B.tank_r, B.ox_tank_h, B.fu_tank_h, B.ox_baffle.Nb, B.fu_baffle.Nb];
    now_  = [VT.radius, VT.ox.height, VT.fu.height, VT.ox.baffle.rings, VT.fu.baffle.rings];
    if any(abs(built - now_) > 1e-9)
        warning('baffle_lut.mat was generated for different tank geometry or baffles; rerun generateBaffleLUT.');
    end
    clear B VT built now_
end

%prop_system_testing();

    