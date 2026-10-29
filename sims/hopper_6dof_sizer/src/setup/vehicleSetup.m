function V = vehicleSetup(IN)
%VEHICLESETUP  The vehicle variables the models read, from the vehicle config.
%   V = vehicleSetup(IN)   IN from mission_inputs (one Monte Carlo
%                          scenario's, or the nominal); its vehicle data
%                          came from vehicles/<vehicle>/<vehicle>_config.m
%
%   Replaces running the sizer during sim setup (src/sizing is now a design
%   tool, sizeVehicle). The structs keep the names and fields the models
%   and setup scripts read, filled from the config:
%     V.TANKS    .singular.radius / .oxidizer / .fuel (.mass, .h) /
%                .total_mass, .volumes
%     V.STRUCT   .legs.k / .legs.mass, .intertank.mass, .total
%     V.ENGINE   .mass, .tvc_natural_freq
%     V.VEH      .mass.dry (x IN.mass_factor), .mass.wet, .k_engine
%     V.OUT      Outputs(...): mass and TWR summary
%     V.r_leg    3 x 4 foot positions; V.r_leg_1..4 rows (older models)
%   New code should read IN.vehicle (the config) instead.

VC = IN.vehicle;

TANKS.singular.radius          = IN.tanks.max_radius;
TANKS.singular.oxidizer.mass   = VC.tanks.ox.mass;
TANKS.singular.oxidizer.h      = IN.tanks.ox_h_as_built;
TANKS.singular.fuel.mass       = VC.tanks.fu.mass;
TANKS.singular.fuel.h          = IN.tanks.fu_h_as_built;
TANKS.singular.total_mass      = VC.tanks.ox.mass + VC.tanks.fu.mass;
TANKS.volumes.oxidizer         = IN.propulsion.ox_vol;
TANKS.volumes.fuel             = IN.propulsion.fu_vol;

STRUCT.legs.k          = VC.legs.stiffness;
STRUCT.legs.mass       = VC.legs.mass;
STRUCT.intertank.mass  = VC.legs.intertank_mass;
STRUCT.total           = STRUCT.legs.mass + STRUCT.intertank.mass;

ENGINE.mass             = VC.engine.mass;
ENGINE.tvc_natural_freq = 1.8 / IN.tvc.rise_time;

% Dry mass from the mass budget (legs and intertank included), dispersed
% by a Monte Carlo scenario's mass_factor
VEH.mass.dry = IN.mass_factor * IN.mass.dry_override;
VEH.mass.wet = VEH.mass.dry + IN.propulsion.oxidizer_mass + IN.propulsion.fuel_mass;
VEH.k_engine = (IN.engine_cg / IN.mission.max_vehicle_height)^2;   % stopgap bending-mode factor (VC.structure)

OUT = Outputs(IN, VEH, TANKS, STRUCT);

r_leg = VC.legs.r_foot;

V.TANKS = TANKS; V.STRUCT = STRUCT; V.ENGINE = ENGINE; V.VEH = VEH; V.OUT = OUT;
V.r_leg = r_leg;
for k = 1:size(r_leg, 2)
    V.(sprintf('r_leg_%d', k)) = r_leg(:, k)';
end
end
