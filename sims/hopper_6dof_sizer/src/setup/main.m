% Vehicle setup (first step of sim_setup): the mission inputs, with the
% vehicle from its config file (vehicles/<vehicle>/<vehicle>_config.m),
% and the vehicle variables the models read. Nothing is sized here; for a
% sizing trade study run sizeVehicle (src/sizing).
IN = mission_inputs(hopperScenario());   % [] = nominal

V = vehicleSetup(IN);
TANKS  = V.TANKS;
STRUCT = V.STRUCT;
ENGINE = V.ENGINE;
VEH    = V.VEH;
OUT    = V.OUT;
r_leg   = V.r_leg;
r_leg_1 = V.r_leg_1;   % per-leg names still read by the v2/v3 models
r_leg_2 = V.r_leg_2;
r_leg_3 = V.r_leg_3;
r_leg_4 = V.r_leg_4;
clear V

fprintf('\n===== VEHICLE: %s (%s) =====\n', IN.vehicle.name, IN.vehicle.revision);
fprintf('Dry Mass: %.2f kg\n', OUT.Vehicle.DryMass);
fprintf('Wet Mass: %.2f kg\n', OUT.Vehicle.WetMass);
fprintf('Propellant: %.2f kg ox, %.2f kg fuel\n', IN.propulsion.oxidizer_mass, IN.propulsion.fuel_mass);
fprintf('Initial TWR: %.2f   Final TWR: %.2f\n', OUT.Vehicle.InitialTWR, OUT.Vehicle.FinalTWR);
fprintf('Thrust: %.1f to %.1f N\n', OUT.Vehicle.minThrust, OUT.Vehicle.maxThrust);
fprintf('Tanks: r %.4f m, ox %.4f m at %.3f m, fuel %.4f m at %.3f m\n', TANKS.singular.radius, ...
    TANKS.singular.oxidizer.h, IN.ox_tank_cg, TANKS.singular.fuel.h, IN.fu_tank_cg);
