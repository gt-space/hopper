function IN = ctv_inputs(IN)
%CTV_INPUTS  Controls test vehicle: overrides on top of the hopper inputs.
%   IN = ctv_inputs(mission_inputs()) is what mission_inputs returns when
%   hopperVehicle is 'ctv' (variantConfig('VEHICLE', 'ctv')).
%
%   The CTV is a ~2.6 kg electric vehicle: one gimbaled EDF (~35 N) for
%   thrust and pitch/yaw, four small props for roll. It flies the same
%   actuator model as the hopper and burns no propellant, so it runs with
%   SLOSH=off and MASS_PROPS=constant (variantConfig's CTV defaults).
%
%   Mass, CG and inertia come from CAD (as in ctv/inputs/mission_inputs.m).
%   Values marked ASSUMPTION are placeholders until measured; check them
%   before trusting a CTV result.

% ---- CAD mass properties (CAD frame: origin at the base centre, Y up) --
IN.ctv.mass = 2614.88 / 1000;                       % kg
IN.ctv.cg_cad = [-41.76; 349.21; 6.07] * 1e-3;      % m, CAD X Y Z
IN.ctv.I_cad = [103179574.33   1374767.18     68051.49;
                  1374767.18  12321923.43  -1150369.65;
                    68051.49  -1150369.65  103816055.97] * 1e-9;  % kg m^2 about the CG, CAD axes
IN.ctv.thrust = 7.875 * 4.44822;                    % N, EDF static thrust

% Body x is the thrust axis (CAD Y); y, z follow right-handed (CAD Z, CAD X)
P = [0 1 0; 0 0 1; 1 0 0];                          % body = P * CAD
cg_body = P * IN.ctv.cg_cad;
IN.mass.dry_cg = cg_body(1);                        % m, axial CG from the base
IN.mass.dry_I = P * IN.ctv.I_cad * P';              % kg m^2, body axes
IN.ctv.cg_lateral = cg_body(2:3);                   % m; the sim uses the axial CG only

% ---- Propellant: none. The propellant model needs tanks that aren't
% empty (the depletion cutoff stops thrust at zero), so it carries 1 g
% placeholders inside the CAD mass; MASS_PROPS=constant keeps them fixed.
IN.propulsion.burns_propellant = false;
IN.propulsion.oxidizer_mass = 1e-3;                 % kg, placeholder
IN.propulsion.fuel_mass = 1e-3;                     % kg, placeholder
IN.mass.dry_override = IN.ctv.mass - IN.propulsion.oxidizer_mass - IN.propulsion.fuel_mass;

% ---- EDF and gimbal ----
IN.propulsion.nominal_thrust = IN.ctv.thrust;      % N
IN.propulsion.throttle_range = [0.2 1.0];          % ASSUMPTION: EDF idle to full
IN.tvc.max_gimbal_angle = deg2rad(8);              % from ctv/inputs/mission_inputs.m
IN.engine_mass = 0.35;                             % kg, ASSUMPTION: EDF + gimbal (inside the CAD mass)
IN.engine_cg = 0.08;                               % m, ASSUMPTION: EDF CG above the base
IN.mount_cg = 0.15;                                % m, ASSUMPTION: gimbal pivot above the base
IN.mission.max_vehicle_height = 0.75;              % m, ASSUMPTION

% ---- Roll props ----
IN.rcs.max_thrust = 16;                            % N, TODO from ctv inputs: replace with datasheet thrust

% ---- Landing legs ----
IN.legs.feet_radial_dist = 0.15;                   % m, ASSUMPTION

% ---- Aerodynamics (Wind_Forces) ----
IN.aero.area = 0.15;                               % m^2, ASSUMPTION: side area
IN.aero.cd = 1.2;                                  % ASSUMPTION
IN.aero.cp = 0.35;                                 % m, ASSUMPTION: centre of pressure from the base

% ---- Wind: calm. ASSUMPTION: the CTV flies in calm air; the hopper's
% 10 m/s default puts ~11 N of side force (~40 % of weight) on it
IN.wind.uwind = 0;                                 % m/s
IN.wind.vwind = 0;                                 % m/s

% ---- Mission: 4 m hop with limits scaled to the CTV ----
IN.mission.target_altitude = 4;                    % m  [req]
IN.mission.max_ascent_vel = 2;                     % m/s [req] ASSUMPTION (CTV inputs had 26)
IN.mission.max_descent_vel = 2;                    % m/s [req] ASSUMPTION
IN.mission.trajectory.apex = 4.5;                  % m
IN.mission.trajectory.t_up = 6;                    % s (peak 1.4 m/s)
IN.mission.trajectory.t_hover = 2;                 % s
IN.mission.trajectory.t_down = 6;                  % s
IN.mission.trajectory.t_term = 2;                  % s
IN.mission.trajectory.z_td = 0.2;                  % m
IN.mission.trajectory.v_td = 0.1;                  % m/s
IN.mission.trajectory.burn_mass = 0;               % kg
end
