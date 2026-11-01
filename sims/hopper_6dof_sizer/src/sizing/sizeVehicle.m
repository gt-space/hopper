function [VC, S] = sizeVehicle(vehicle, varargin)
%SIZEVEHICLE  Conceptual sizing trade study; can draft a vehicle config.
%   sizeVehicle()                    size the current vehicle (hopperVehicle)
%                                    from its mission inputs and print the
%                                    sizer's numbers next to the config's
%   [VC, S] = sizeVehicle('hopper')  VC: the vehicle's config with the sized
%                                    values in place; S: the sizer outputs
%   sizeVehicle(vehicle, 'Write', file)
%                                    also write VC as a draft config file
%                                    (writeVehicleConfig) for review
%
%   The sizer (size_tanks, size_structures, size_engine, size_avionics) no
%   longer runs in sim setup: the sim flies the config file. Use this for
%   design trades, then copy the values you adopt into
%   vehicles/<vehicle>/<vehicle>_config.m (or start from the draft).
%   It sizes: empty tank masses (pressure-vessel walls), leg stiffness and
%   mass (buckling-sized rods), intertank mass. CGs, sensors and anything
%   from CAD come from the current config unchanged.

if nargin < 1 || isempty(vehicle), vehicle = hopperVehicle(); end
vehicle = char(vehicle);
p = inputParser;
p.addParameter('Write', '', @(x) ischar(x) || isstring(x));
p.parse(varargin{:});

old = hopperVehicle(vehicle);
restore = onCleanup(@() hopperVehicle(old));
IN = mission_inputs();

S.TANKS  = size_tanks(IN);
S.AVI    = size_avionics(IN);
S.ENGINE = size_engine(IN);
VEH.mass.dry = IN.mass.dry_override;
VEH.mass.wet = VEH.mass.dry + IN.propulsion.oxidizer_mass + IN.propulsion.fuel_mass;
[S.STRUCT, ~] = size_structures(IN, VEH);

VC = feval([vehicle '_config']);
cur = VC;
VC.tanks.ox.mass = S.TANKS.singular.oxidizer.mass;
VC.tanks.fu.mass = S.TANKS.singular.fuel.mass;
VC.legs.stiffness = S.STRUCT.legs.k;
VC.legs.mass = S.STRUCT.legs.mass;
VC.legs.intertank_mass = S.STRUCT.intertank.mass;
VC.revision = sprintf('sizer draft %s', char(datetime('today', 'Format', 'yyyy-MM-dd')));

fprintf('\n===== SIZING: %s (mission inputs + config) =====\n', vehicle);
row = @(name, unit, sized, now_) fprintf('%-26s %12.4f %12.4f  %s\n', name, sized, now_, unit);
fprintf('%-26s %12s %12s\n', '', 'sizer', 'config');
row('Ox tank mass', 'kg', VC.tanks.ox.mass, cur.tanks.ox.mass);
row('Fuel tank mass', 'kg', VC.tanks.fu.mass, cur.tanks.fu.mass);
row('Ox tank height (V/pi r^2)', 'm', S.TANKS.singular.oxidizer.h, cur.tanks.ox.height);
row('Fuel tank height', 'm', S.TANKS.singular.fuel.h, cur.tanks.fu.height);
row('Leg stiffness', 'MN/m', VC.legs.stiffness / 1e6, cur.legs.stiffness / 1e6);
row('Legs mass', 'kg', VC.legs.mass, cur.legs.mass);
row('Engine mass', 'kg', S.ENGINE.mass, cur.engine.mass);
row('Avionics mass', 'kg', S.AVI.mass, cur.components.avionics.mass);
fprintf('\nTanks, clustered (3): %.2f kg   concentric: %.2f kg   stacked: %.2f kg\n', ...
    S.TANKS.clustered.total_mass, S.TANKS.concentric.total_mass, S.TANKS.singular.total_mass);
fprintf('Engine: injector %.2f kg, TCA %.2f kg, TVC %.2f kg\n', S.ENGINE.inj_mass, S.ENGINE.TCA_mass, S.ENGINE.TVC_mass);
fprintf('Battery: %.2f Ah, %d cells\n', S.AVI.battery_capacity, S.AVI.num_cells);
fprintf('Heights are not drafted: the config keeps the as-built tank lengths.\n');

if strlength(p.Results.Write) > 0
    sized = {'tanks.ox.mass', 'tanks.fu.mass', 'legs.stiffness', 'legs.mass', 'legs.intertank_mass'};
    writeVehicleConfig(VC, p.Results.Write, 'Sized', sized);
end
end
