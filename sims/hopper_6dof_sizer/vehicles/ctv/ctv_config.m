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
% at that file's placeholder thrust. TODO: datasheet value (typical 65 mm
% props give ~0.4 N). Roll torque limit: max_thrust x arm.
VC.rcs.max_thrust = 16;                        % N, placeholder
VC.rcs.arm        = 4.25*0.0254;               % m, ASSUMPTION: still the hopper's tank radius; measure the prop arm
% The hopper's 0.1 N m quantization is several times the torque the CTV's
% roll loop commands (a few mN m), so the props never moved
VC.rcs.resolution = 1e-3;                      % N m, ASSUMPTION: ~0.01 N of prop thrust

% ---- Random roll disturbance torque (the harness's Random Gusts) ----
% The hopper's, scaled by roll inertia squared: the same angular
% acceleration (20 (N m)^2 x (0.012329 / 12.6948)^2). The hopper's own
% sigma, 4.5 N m, is 2.6x the CTV's whole roll authority. The CTV's real
% roll disturbance is the EDF's reaction torque, which is not modelled.
VC.disturbance.roll_torque_var = 1.886e-5;     % (N m)^2, ASSUMPTION

% ---- Landing legs ----
VC.legs.foot_radius = 0.15;                    % m, ASSUMPTION
VC.legs.stiffness   = 2186514.3576498032;      % N/m per leg: sizer 2026-10, at the CTV's mass
VC.legs.mass        = 0.44908343062561118;     % kg, all legs: sizer 2026-10

% ---- Structure: no bending-mode coupling ----
% The hopper's stopgap bending mode (k_engine = (engine cg_x / height)^2,
% at 10 Hz, inherited from hopper_config) was sized for a 2.44 m vehicle.
% The CTV is a stiff 0.75 m frame with no first mode near its control
% bandwidth, so it flies rigid until a structural model says otherwise.
VC.structure.modal.enabled = false;

% ---- Aerodynamics (Wind_Forces) ----
VC.aero.area = 0.15;                           % m^2, ASSUMPTION: side area
VC.aero.cd   = 1.2;                            % ASSUMPTION
VC.aero.cp_x = 0.35;                           % m, ASSUMPTION: centre of pressure

% ---- Navigation filter tuning ----
% The hopper's torque PSD scaled by inertia (same angular-acceleration
% noise per axis, J roll 0.0123, pitch/yaw 0.105 kg m^2): the hopper's
% 1 (N m)^2/Hz on this inertia lets the roll rate and gyro bias trade
% freely and the attitude diverges within 1.5 s
VC.nav.torque_psd = [1.4e-4 3e-6 3e-6];        % (N m)^2/Hz, roll pitch yaw

% ---- Sensors: the hopper's lidar layout; GPS antennas near the top ----
% Antenna 1 where CTV_GMEKF's gmekf_correct assumes it: [0.1; 0; -0.2] m from
% the CG in its frame, i.e. 0.2 m above the CG and 0.1 m along body z
% (the CTV keeps its CG); antenna 2 is unused (no differential GPS)
VC.sensors.gps.antenna = [VC.mass.dry_cg + 0.2, VC.mass.dry_cg + 0.2; 0 0; 0.1 -0.1];   % m, ASSUMPTION: body x y z from the base
end
