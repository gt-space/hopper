%% Hopper Vehicle Sizer - Main Entry Point
clear; clc;

%% =======================
% Add project paths
%% =======================
root = fileparts(mfilename('fullpath'));
addpath(genpath(root));

%% =======================
% Load inputs
%% =======================
IN = mission_inputs();

%% =======================
% User-defined fixed inputs
%% =======================
IN.propulsion.Isp_vac = 245;              % s (guess)
IN.propulsion.OF = 4.0;
IN.propulsion.chamber_pressure = 24e5;    % Pa (300 psi)

IN.avionics.mass = 3.0;                   % kg
IN.avionics.power = 80;                   % W

%% =======================
% Initial vehicle mass guess
%% =======================
vehicle = struct();
vehicle.mass_dry = 40;    % kg initial guess
vehicle.mass_wet = vehicle.mass_dry;

%% =======================
% Convergence settings
%% =======================
tol = 0.5;          % kg
max_iter = 20;

fprintf('Starting vehicle sizing loop...\n');

for iter = 1:max_iter

    fprintf('\nIteration %d\n', iter);

    %% =======================
    % Engine sizing
    %% =======================
    engine = size_engine(IN, vehicle);

    %% =======================
    % Propellant sizing
    %% =======================
    prop = size_propellant(IN, engine, vehicle);

    %% =======================
    % Tank sizing
    %% =======================
    tanks = size_tanks(IN, prop);

    %% =======================
    % Update vehicle masses
    %% =======================
    m_prop = prop.m_total;
    m_tanks = tanks.ox.mass + tanks.fu.mass;

    vehicle.mass_dry = ...
        IN.structures.payload_mass + ...
        IN.avionics.mass + ...
        engine.mass + ...
        m_tanks;

    vehicle.mass_wet = vehicle.mass_dry + m_prop + IN.margins.mass_growth;

    %% =======================
    % Structural sizing
    %% =======================
    structures = size_structures(IN, vehicle);
    vehicle.mass_dry = vehicle.mass_dry + structures.total;
    vehicle.mass_wet = vehicle.mass_wet + structures.total;

    %% =======================
    % Battery sizing
    %% =======================
    battery = size_batteries(IN, []);
    vehicle.mass_dry = vehicle.mass_dry + battery.mass;
    vehicle.mass_wet = vehicle.mass_wet + battery.mass;

    %% =======================
    % Convergence check
    %% =======================
    if iter > 1
        delta = abs(vehicle.mass_wet - mass_prev);
        fprintf('  Wet mass = %.2f kg (Δ = %.2f)\n', vehicle.mass_wet, delta);
        if delta < tol
            fprintf('Converged.\n');
            break;
        end
    else
        fprintf('  Wet mass = %.2f kg\n', vehicle.mass_wet);
    end

    mass_prev = vehicle.mass_wet;

end

%% =======================
% Final outputs
%% =======================
Vehicle = vehicle;

OUT = struct();
OUT.vehicle.mass_dry = vehicle.mass_dry;
OUT.vehicle.mass_wet = vehicle.mass_wet;
OUT.vehicle.TWR_initial = engine.thrust_nom / (vehicle.mass_wet * IN.const.g0);
OUT.vehicle.TWR_final = engine.thrust_nom / (vehicle.mass_dry * IN.const.g0);

OUT.propellant = prop;
OUT.engine = engine;
OUT.tanks = tanks;
OUT.structures = structures;
OUT.battery = battery;

%% =======================
% Print summary
%% =======================
fprintf('\n=== FINAL VEHICLE SUMMARY ===\n');
fprintf('Dry Mass: %.2f kg\n', vehicle.mass_dry);
fprintf('Wet Mass: %.2f kg\n', vehicle.mass_wet);
fprintf('Initial TWR: %.2f\n', OUT.vehicle.TWR_initial);
fprintf('Final TWR: %.2f\n', OUT.vehicle.TWR_final);
%% =======================
% Save results
%% =======================
save('hopper_sizing_output.mat','IN','OUT','Vehicle');

