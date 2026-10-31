function VC = hopper_config()
%HOPPER_CONFIG  The hopper as designed: geometry, masses and hardware.
%   VC = hopper_config()   read through vehicleConfig('hopper'), which
%                          checks it and adds the derived values
%
%   This file is the source of truth for the vehicle. Change a value here
%   and every subsystem sees it (docs/VEHICLE_CONFIG.md). Requirements,
%   the trajectory, wind and Monte Carlo dispersions are not vehicle data;
%   they stay in inputs/mission_inputs.m.
%
%   Frame: body x up the thrust axis from the base reference point (the
%   bottom centre of the vehicle, x = 0 on the pad); y and z complete a
%   right-handed set. Every position below is an x in this frame (m).
%
%   Values tagged "sizer 2026-10" are the conceptual sizer's numbers
%   (src/sizing), frozen when the sizer stopped running in sim setup.
%   Replace them with as-built values.

VC.name     = 'hopper';
VC.revision = 'CAD 2026-10-29';
VC.source   = ['dry block: SOLIDWORKS mass properties of Updated Vehicle Assy (2026-10-29); ' ...
               'other CAD masses and CGs (mission_inputs.m 2026-10); sizer values where tagged'];
VC.height   = 2.44;                            % m, base to top

% ---- Dry mass: CAD mass properties ----
% "Updated Vehicle Assy - Mass Properties" (SOLIDWORKS IMassProperty2,
% 2026-10-29): 552 components, hidden ones included. Replaces the 116.86 kg
% mass-budget stand-in, whose CG and inertia came from the component list.
% CAD frame: assembly origin, Y up the thrust axis. The footpad bottoms (the
% assembly's lowest point, swBoundingBox) are at CAD Y = -1.208957 m; that
% is body x = 0.
% In the CAD: legs, intertank structure, thrust structure, gimbal, TVC
% actuators, engine, tanks and baffles, COPVs with their gas, fluids
% hardware, the avionics enclosure. Not in it: propellant (added in flight),
% avionics boards and battery, the payload, the EDF fans. The components
% block below is not added to this (mass_properties_setup uses dry_cg/dry_I).
VC.cad.mass   = 114.894;                       % kg
VC.cad.cg     = [-0.003159; -0.169015; 0.003797];   % m, CAD X Y Z from the assembly origin
VC.cad.I      = [51.5454   0.2508  -0.0406;
                  0.2508  12.6948   0.1037;
                 -0.0406   0.1037  52.3340];  % kg m^2 about the CG, CAD axes, as exported
y_base = -1.208957;                            % m, CAD Y of the footpad bottoms
% Body x is the thrust axis (CAD Y); y, z follow right-handed (CAD Z, CAD X)
VC.cad.to_body = [0 1 0; 0 0 1; 1 0 0];       % body = to_body * CAD
cg_body = VC.cad.to_body * (VC.cad.cg - [0; y_base; 0]);
I_body  = VC.cad.to_body * VC.cad.I * VC.cad.to_body';

VC.mass.dry = VC.cad.mass;                     % kg
VC.mass.dry_cg = cg_body(1);                   % m, axial CG from the base (1.040)
VC.mass.cg_lateral = cg_body(2:3);             % m (3.8 mm, -3.2 mm); the sim uses the axial CG only
VC.mass.I_full = I_body;                       % kg m^2, body axes, full CAD tensor
% Principal moments only, as for the CTV: the products are < 0.5 % of the
% pitch/yaw inertia, and the principal axes are within 3 deg of the body axes.
VC.mass.dry_I = diag(diag(I_body));            % kg m^2, [roll pitch yaw] = [12.69 52.33 51.55]

% ---- Propellant tanks (stacked, LOx above RP-1) ----
VC.tanks.geometry = 'stacked';
VC.tanks.material = 'Al6061';
VC.tanks.radius   = 4.25*0.0254;               % m, inner radius (4.25 in)
VC.tanks.ox.height    = 29.18*0.0254;          % m, as-built length (29.18 in)
VC.tanks.ox.cg_x      = 1.854;                 % m, tank centre
VC.tanks.ox.wall      = 0.003175;              % m
VC.tanks.ox.mass      = 4.6940418420227319;    % kg, empty tank: sizer 2026-10
VC.tanks.ox.volume    = 0.01638;               % m^3
VC.tanks.ox.prop_mass = 16;                    % kg, LOx load
VC.tanks.ox.density   = 1141;                  % kg/m^3
VC.tanks.ox.viscosity = 0.17e-6;               % m^2/s, kinematic
% Ring baffles (generateBaffleLUT; sweepBaffleDesign 2026-10-20), over the full tank height
VC.tanks.ox.baffle    = struct('rings', 15, 'width', 0.5 * VC.tanks.radius, 'thickness', 0.5e-3);   % -, m, m
VC.tanks.fu.height    = 18.91*0.0254;          % m, as-built length (18.91 in)
VC.tanks.fu.cg_x      = 1.0668;                % m, tank centre
VC.tanks.fu.wall      = 0.003175;              % m
VC.tanks.fu.mass      = 5.6082853202836009;    % kg, empty tank: sizer 2026-10
VC.tanks.fu.volume    = 0.02034;               % m^3
VC.tanks.fu.prop_mass = 13;                    % kg, RP-1 load
VC.tanks.fu.density   = 800;                   % kg/m^3
VC.tanks.fu.viscosity = 1.3e-6;                % m^2/s, kinematic
VC.tanks.fu.baffle    = struct('rings', 11, 'width', 0.5 * VC.tanks.radius, 'thickness', 0.5e-3);   % -, m, m

% ---- Engine and TVC ----
VC.engine.mass           = 37.83;              % kg, CAD (engine, TVC, feed hardware)
VC.engine.cg_x           = 0.3556;             % m
VC.engine.gimbal_x       = 0.45;               % m, gimbal pivot
VC.engine.offset         = [0 0];              % m, thrust-line offset [y z] (Monte Carlo disperses it)
VC.engine.nominal_thrust = 500 * 4.44822;      % N (500 lbf)
VC.engine.throttle_range = [0.4 1.1];          % fraction of nominal
VC.engine.throttle.rate_limit = 180 * 4.44822; % N/s, throttle valve
VC.engine.throttle.latency    = 0.01;          % s
VC.engine.throttle.resolution = 5;             % N
VC.engine.tvc.max_angle  = deg2rad(15);        % rad
VC.engine.tvc.max_rate   = deg2rad(15);        % rad/s
VC.engine.tvc.latency    = 0.01;               % s
VC.engine.tvc.rise_time  = 0.05;               % s, actuator 10-90 % rise
VC.engine.tvc.damping    = 0.5;                % actuator damping ratio
VC.engine.tvc.resolution = deg2rad(0.1);       % rad

% ---- Roll control (EDF fans) ----
VC.rcs.max_thrust = 16;                        % N per fan (roll torque limit: max_thrust x arm)
VC.rcs.arm        = VC.tanks.radius;           % m, fan moment arm (fans at the tank wall)
VC.rcs.tau        = 0.1;                       % s, fan spin-up time constant
VC.rcs.resolution = 0.1;                       % N m, command quantization
VC.rcs.throttle_rate = 50;                     % N/s

% ---- Random roll disturbance torque (the harness's Random Gusts) ----
% On 5 % of 0.1 s intervals, a Gaussian roll torque redrawn every second
VC.disturbance.roll_torque_var = 20;           % (N m)^2

% ---- Landing legs ----
VC.legs.n           = 4;
VC.legs.foot_radius = 0.5;                     % m, feet from the axis
VC.legs.azimuth     = [0 180 90 270];         % deg about body x from +y, in leg order (legs 1-4: +y -y +z -z)
VC.legs.stiffness   = 16970822.627420306;      % N/m per leg: sizer 2026-10 (buckling-sized rod)
VC.legs.damping     = 1.0;                     % strut damping constant (x sqrt(wet mass))
VC.legs.mass        = 3.4856003663531885;      % kg, all legs: sizer 2026-10
VC.legs.intertank_mass = 3;                    % kg, sizer 2026-10 (guess)

% ---- Other components (the pre-CAD dry block; not used while dry_cg is set) ----
VC.components.structures = struct('mass', 34.26, 'cg_x', 1.193);    % kg, m (CAD)
VC.components.avionics   = struct('mass', 11.52, 'cg_x', 1.4224);   % kg, m (CAD)
VC.components.fluids     = struct('mass', 28.41, 'cg_x', 1.0668);   % kg, m (CAD)
VC.components.payload    = struct('mass', 15,    'cg_x', 2);        % kg, m (CPLC requirement)
% COPVs: mass and gas mass as the mass model sums them, one bottle's geometry
VC.components.copv = struct('mass', 11, 'gas_mass', 6, 'amount', 3, ...
    'r', 0.18, 't', 13.21/1000, 'h', 0.56, 'cg_x', 1.21, 'cg_y', 0.23);  % kg, kg, -, m, m, m, m, m

% ---- Structure (placeholder for a structural response model) ----
% The first bending mode feeds the TVC angle through k_engine = (engine
% cg_x / height)^2 (Actuators_TVC_RCS). A stopgap, not a design value.
VC.structure.modal.natural_frequency = 10;     % Hz
VC.structure.modal.damping           = 0.01;

% ---- Aerodynamics (Wind_Forces) ----
VC.aero.area = 0.5216;                         % m^2, side area
VC.aero.cd   = 8.494;
VC.aero.cp_x = 0.9;                            % m, centre of pressure

% ---- Navigation filter tuning that scales with the vehicle ----
% Torque disturbance PSD per body axis [roll pitch yaw] (the GMEKF's Qc
% torque block): the commanded torque leaves out wind, slosh and other
% moments, ~7 N m here. Scales with inertia: keep (torque / J)^2 when
% porting to another vehicle (see ctv_config).
VC.nav.torque_psd = [1 1 1];                   % (N m)^2/Hz

% ---- Sensor mounting (sim_sensor_params) ----
% navigation: four lidars 0.25 m off-axis at 0/90/180/270 deg, canted 20 deg
VC.sensors.lidar.azimuth = [0 90 180 270];     % deg, about body x
VC.sensors.lidar.cant    = 20;                 % deg from the long axis
VC.sensors.lidar.radius  = 0.25;               % m off-axis
VC.sensors.lidar.mount_x = 0.1;                % m (ASSUMPTION)
% GPS antennas (dual, for differential GPS heading), columns 1 and 2, from
% the base. ASSUMPTION: the nav team's placeholder, 0.25 m either side of
% the axis and 0.2 m above the nominal CG; replace with CAD
VC.sensors.gps.antenna = [1.343 1.343; 0 0; 0.25 -0.25];   % m, body x y z
VC.sensors.imu.period  = 1/500;                % s, ADIS16500 at the filter rate (nav team: DEC_RATE = 3)
end
