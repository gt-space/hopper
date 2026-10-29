function VC = ctv_config()
%CTV_CONFIG  The controls test vehicle: geometry, masses and hardware.
%   VC = ctv_config()   read through vehicleConfig('ctv')
%
%   The CTV is a ~2.6 kg electric vehicle: one gimbaled EDF (~35 N) for
%   thrust and pitch/yaw, four small props for roll. It burns no
%   propellant, so it runs with SLOSH=off and MASS_PROPS=constant
%   (variantConfig's CTV defaults). Same frame as hopper_config: body x up
%   the thrust axis from the bottom centre.
%
%   It starts from the hopper's file, so any field not set here is still a
%   hopper value: the tanks and the component list. Neither moves the CTV
%   (its mass properties come from CAD and slosh is off), but they are not
%   CTV data. Values marked ASSUMPTION are placeholders until measured;
%   check them before trusting a CTV result.

VC = hopper_config();
VC.name     = 'ctv';
VC.revision = 'CAD 2026-10';
VC.source   = 'CAD mass properties (ctv/inputs/mission_inputs.m); ASSUMPTION where tagged';
VC.height   = 0.75;                            % m, ASSUMPTION

% ---- CAD mass properties (CAD frame: origin at the base centre, Y up) ----
VC.cad.mass   = 2614.88 / 1000;                % kg
VC.cad.cg     = [-41.76; 349.21; 6.07] * 1e-3; % m, CAD X Y Z
VC.cad.I      = [103179574.33   1374767.18     68051.49;
                   1374767.18  12321923.43  -1150369.65;
                     68051.49  -1150369.65  103816055.97] * 1e-9;  % kg m^2 about the CG, CAD axes
% Body x is the thrust axis (CAD Y); y, z follow right-handed (CAD Z, CAD X)
VC.cad.to_body = [0 1 0; 0 0 1; 1 0 0];       % body = to_body * CAD
cg_body = VC.cad.to_body * VC.cad.cg;
I_body  = VC.cad.to_body * VC.cad.I * VC.cad.to_body';

% ---- Propellant: none. The propellant model needs tanks that aren't
% empty (the depletion cutoff stops thrust at zero), so it carries 1 g
% placeholders inside the CAD mass; MASS_PROPS=constant keeps them fixed.
VC.burns_propellant = false;
VC.tanks.ox.prop_mass = 1e-3;                  % kg, placeholder
VC.tanks.fu.prop_mass = 1e-3;                  % kg, placeholder

VC.mass.dry = VC.cad.mass - VC.tanks.ox.prop_mass - VC.tanks.fu.prop_mass;   % kg
VC.mass.dry_cg = cg_body(1);                   % m, axial CG from the base
VC.mass.cg_lateral = cg_body(2:3);             % m; the sim uses the axial CG only
VC.mass.I_full = I_body;                       % kg m^2, body axes, full CAD tensor
% Principal moments only, as with the CG: the CAD products of inertia
% (up to 9 % of the small roll inertia) make the LQR's flight drift 8.5 m
% and tilt 20 deg in calm air, with no disturbance at all. Investigate that
% coupling before trusting a CTV flight with the full tensor.
VC.mass.dry_I = diag(diag(I_body));            % kg m^2

% ---- EDF and gimbal ----
VC.engine.nominal_thrust = 7.875 * 4.44822;    % N, EDF static thrust
VC.engine.throttle_range = [0.2 1.0];          % ASSUMPTION: EDF idle to full
VC.engine.mass     = 0.35;                     % kg, ASSUMPTION: EDF + gimbal (inside the CAD mass)
VC.engine.cg_x     = 0.08;                     % m, ASSUMPTION: EDF CG above the base
VC.engine.gimbal_x = 0.15;                     % m, ASSUMPTION: gimbal pivot above the base
VC.engine.tvc.max_angle = deg2rad(8);          % rad, from ctv/inputs/mission_inputs.m
VC.engine.tvc.max_rate  = deg2rad(200);        % rad/s, ASSUMPTION: hobby servo (hopper 15 deg/s)
% EDF throttle: an ESC, not the hopper's throttle valve (5 N steps there
% are 20 % of the CTV's weight and made thrust bang-bang)
VC.engine.throttle.resolution = 0.05;          % N, ASSUMPTION

% ---- Roll props ----
% Four 65 mm props on coreless motors (ctv/inputs/mission_inputs.m), still
% at that file's placeholder thrust. TODO: datasheet value. The roll
% authority matters: with the current LQR the CTV holds attitude under a
% 1 mm engine offset or 1.5 m/s wind at 16 N but tumbles at 2 N or less
% (typical 65 mm props give ~0.4 N), so check this before flying.
VC.rcs.max_thrust = 16;                        % N, placeholder
VC.rcs.arm        = 4.25*0.0254;               % m, ASSUMPTION: still the hopper's tank radius; measure the prop arm

% ---- Landing legs ----
VC.legs.foot_radius = 0.15;                    % m, ASSUMPTION
VC.legs.stiffness   = 2186514.3576498032;      % N/m per leg: sizer 2026-10, at the CTV's mass
VC.legs.mass        = 0.44908343062561118;     % kg, all legs: sizer 2026-10

% ---- Aerodynamics (Wind_Forces) ----
VC.aero.area = 0.15;                           % m^2, ASSUMPTION: side area
VC.aero.cd   = 1.2;                            % ASSUMPTION
VC.aero.cp_x = 0.35;                           % m, ASSUMPTION: centre of pressure
end
