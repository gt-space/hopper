function IN = vehicleInputs(IN, VC)
%VEHICLEINPUTS  Put a vehicle config's values where the setup code reads them.
%   IN = vehicleInputs(IN, VC)   called by mission_inputs, before the
%                                Monte Carlo dispersions
%
%   The setup scripts and models read vehicle data from IN fields that
%   predate the config file (IN.engine_cg, IN.tanks.max_radius, ...). This
%   fills them from VC, so VC (vehicles/<vehicle>/<vehicle>_config.m) is
%   the only place a vehicle value is written. IN.vehicle is VC itself,
%   for new code: read IN.vehicle.* rather than adding IN fields here.
%   A Monte Carlo dispersion changes the IN field, not IN.vehicle, so
%   IN.vehicle is always the nominal design.

IN.vehicle = VC;

% ---- Mass ----
IN.mission.max_vehicle_height = VC.height;
IN.mass.dry_override = VC.mass.dry;
if isfield(VC.mass, 'dry_cg'), IN.mass.dry_cg = VC.mass.dry_cg; end
if isfield(VC.mass, 'dry_I'),  IN.mass.dry_I  = VC.mass.dry_I;  end
C = VC.components;
IN.structures.payload_mass    = C.payload.mass;
IN.structures.structures_mass = C.structures.mass;
IN.avionics.avi_mass          = C.avionics.mass;
IN.fluids_mass                = C.fluids.mass;
IN.engine_mass                = VC.engine.mass;
IN.press.copv.mass     = C.copv.mass;
IN.press.copv.gas_mass = C.copv.gas_mass;
IN.press.copv.amount   = C.copv.amount;
IN.press.copv.r        = C.copv.r;
IN.press.copv.t        = C.copv.t;
IN.press.copv.h        = C.copv.h;

% ---- CG positions ----
IN.engine_cg          = VC.engine.cg_x;
IN.mount_cg           = VC.engine.gimbal_x;
IN.ox_tank_cg         = VC.tanks.ox.cg_x;
IN.ox_tank_wall_thick = VC.tanks.ox.wall;
IN.fu_tank_cg         = VC.tanks.fu.cg_x;
IN.fu_tank_wall_thick = VC.tanks.fu.wall;
IN.structures_cg      = C.structures.cg_x;
IN.avi_cg             = C.avionics.cg_x;
IN.fluids_cg          = C.fluids.cg_x;
IN.payload_cg         = C.payload.cg_x;
IN.copv_r   = C.copv.r;
IN.copv_t   = C.copv.t;
IN.copv_h   = C.copv.h;
IN.copv_zcg = C.copv.cg_x;
IN.copv_ycg = C.copv.cg_y;

% ---- Propulsion ----
IN.propulsion.burns_propellant = VC.burns_propellant;
IN.propulsion.oxidizer_mass  = VC.tanks.ox.prop_mass;
IN.propulsion.fuel_mass      = VC.tanks.fu.prop_mass;
IN.propulsion.ox_vol         = VC.tanks.ox.volume;
IN.propulsion.fu_vol         = VC.tanks.fu.volume;
IN.propulsion.nominal_thrust = VC.engine.nominal_thrust;
IN.propulsion.throttle_range = VC.engine.throttle_range;
IN.off_axis_y = VC.engine.offset(1);   % thrust-line offset, body y
IN.off_axis_x = VC.engine.offset(2);   % ... body z (the name predates the frame)
IN.throttle_valve.rate_limit = VC.engine.throttle.rate_limit;
IN.throttle_valve_latency    = VC.engine.throttle.latency;
IN.throttle_valve.resolution = VC.engine.throttle.resolution;

% ---- Tanks ----
IN.tanks.geometry      = VC.tanks.geometry;
IN.tanks.material      = VC.tanks.material;
IN.tanks.max_radius    = VC.tanks.radius;
IN.tanks.ox_h_as_built = VC.tanks.ox.height;
IN.tanks.fu_h_as_built = VC.tanks.fu.height;

% ---- TVC and roll control ----
T = VC.engine.tvc;
IN.tvc.max_gimbal_angle = T.max_angle;
IN.tvc.max_gimbal_rate  = T.max_rate;
IN.tvc.actuator_latency = T.latency;
IN.tvc.rise_time        = T.rise_time;
IN.tvc.resolution       = T.resolution;
IN.tvc.damping_ratio    = T.damping;
IN.rcs.max_thrust    = VC.rcs.max_thrust;
IN.rcs.throttle_rate = VC.rcs.throttle_rate;
IN.rcs.resolution    = VC.rcs.resolution;
IN.rcs_time_constant = VC.rcs.tau;

% ---- Legs, structure, aero, sensors ----
IN.legs.feet_radial_dist = VC.legs.foot_radius;
IN.legs.damping_const    = VC.legs.damping;
IN.vehicle_natural_frequency = VC.structure.modal.natural_frequency;
IN.vehicle_damping_ratio     = VC.structure.modal.damping;
IN.aero.area = VC.aero.area;
IN.aero.cd   = VC.aero.cd;
IN.aero.cp   = VC.aero.cp_x;
if isfield(VC, 'sensors')
    IN.sensors = VC.sensors;   % sim_sensor_params merges it over its defaults
end
end
