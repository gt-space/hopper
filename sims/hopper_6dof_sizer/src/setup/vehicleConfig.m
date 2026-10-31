function VC = vehicleConfig(vehicle)
%VEHICLECONFIG  A vehicle's config file, checked, with its derived values.
%   VC = vehicleConfig()          the current vehicle (hopperVehicle)
%   VC = vehicleConfig('ctv')     vehicles/ctv/ctv_config.m
%
%   The config file (vehicles/<vehicle>/<vehicle>_config.m) holds the
%   design values; this adds what follows from them, in one place, so no
%   consumer recomputes them:
%     VC.tanks.ox/fu.floor_x   tank bottom (cg_x - height/2)
%     VC.mass.wet              dry mass + propellant loads
%     VC.legs.r_foot           3 x n foot positions, body frame (legs.azimuth order)
%     VC.engine.tvc.natural_freq   1.8 / rise_time (rad/s)
%     VC.structure.modal.k_engine  (engine cg_x / height)^2, the stopgap
%                                  bending-mode factor on the TVC angle
%   and checks it: required fields, positive masses and sizes, positions
%   inside the vehicle, tanks that don't overlap. See docs/VEHICLE_CONFIG.md.

if nargin < 1 || isempty(vehicle), vehicle = hopperVehicle(); end
vehicle = char(vehicle);
fn = [vehicle '_config'];
if exist(fn, 'file') ~= 2
    error('vehicleConfig:missing', 'No config for vehicle ''%s'' (expected vehicles/%s/%s.m).', vehicle, vehicle, fn);
end
VC = feval(fn);
if ~isfield(VC, 'burns_propellant'), VC.burns_propellant = true; end

check(VC, vehicle);

% ---- Derived values ----
for t = ["ox", "fu"]
    VC.tanks.(t).floor_x = VC.tanks.(t).cg_x - VC.tanks.(t).height / 2;
end
VC.mass.wet = VC.mass.dry + VC.tanks.ox.prop_mass + VC.tanks.fu.prop_mass;
az = deg2rad(VC.legs.azimuth);
VC.legs.r_foot = [zeros(1, VC.legs.n); VC.legs.foot_radius * round(cos(az), 12); ...
                  VC.legs.foot_radius * round(sin(az), 12)];
VC.engine.tvc.natural_freq = 1.8 / VC.engine.tvc.rise_time;
VC.structure.modal.k_engine = (VC.engine.cg_x / VC.height)^2;
end

function check(VC, vehicle)
need = {'name', 'height', 'mass.dry', 'tanks.radius', 'tanks.ox.height', 'tanks.ox.cg_x', ...
        'tanks.ox.mass', 'tanks.ox.prop_mass', 'tanks.fu.height', 'tanks.fu.cg_x', 'tanks.fu.mass', ...
        'tanks.fu.prop_mass', 'engine.mass', 'engine.cg_x', 'engine.gimbal_x', 'engine.nominal_thrust', ...
        'engine.throttle_range', 'engine.tvc.rise_time', 'rcs.max_thrust', 'rcs.arm', 'legs.n', ...
        'legs.foot_radius', 'legs.azimuth', 'legs.stiffness', 'aero.area', 'aero.cp_x', 'structure.modal.natural_frequency', ...
        'disturbance.roll_torque_var'};
bad = {};
for k = 1:numel(need)
    if ~hasPath(VC, need{k}), bad{end+1} = need{k}; end %#ok<AGROW>
end
if ~isempty(bad)
    error('vehicleConfig:field', '%s config is missing %s.', vehicle, strjoin(bad, ', '));
end
pos = {'height', VC.height; 'mass.dry', VC.mass.dry; 'tanks.radius', VC.tanks.radius; ...
       'tanks.ox.height', VC.tanks.ox.height; 'tanks.fu.height', VC.tanks.fu.height; ...
       'engine.mass', VC.engine.mass; 'engine.nominal_thrust', VC.engine.nominal_thrust; ...
       'legs.stiffness', VC.legs.stiffness; 'engine.tvc.rise_time', VC.engine.tvc.rise_time};
for k = 1:size(pos, 1)
    if ~(isnumeric(pos{k, 2}) && isscalar(pos{k, 2}) && pos{k, 2} > 0)
        error('vehicleConfig:value', '%s config: %s must be a positive number.', vehicle, pos{k, 1});
    end
end
tr = VC.engine.throttle_range;
if numel(tr) ~= 2 || tr(1) <= 0 || tr(1) >= tr(2)
    error('vehicleConfig:value', '%s config: engine.throttle_range must be [min max] with 0 < min < max.', vehicle);
end
inside = {'engine.cg_x', VC.engine.cg_x; 'engine.gimbal_x', VC.engine.gimbal_x; 'aero.cp_x', VC.aero.cp_x};
if VC.burns_propellant
    % the tanks only have to fit a vehicle that carries propellant (the
    % CTV keeps the hopper's tanks as unused placeholders)
    inside = [inside; {'ox tank top', VC.tanks.ox.cg_x + VC.tanks.ox.height / 2; ...
                       'fu tank bottom', VC.tanks.fu.cg_x - VC.tanks.fu.height / 2}];
    lo = sort([VC.tanks.ox.cg_x - VC.tanks.ox.height / 2, VC.tanks.fu.cg_x - VC.tanks.fu.height / 2]);
    hi = sort([VC.tanks.ox.cg_x + VC.tanks.ox.height / 2, VC.tanks.fu.cg_x + VC.tanks.fu.height / 2]);
    if hi(1) > lo(2) + 1e-9
        error('vehicleConfig:geometry', '%s config: the ox and fuel tanks overlap.', vehicle);
    end
end
for k = 1:size(inside, 1)
    if inside{k, 2} < 0 || inside{k, 2} > VC.height
        error('vehicleConfig:geometry', '%s config: %s (%.3f m) is outside the vehicle (0 to %.3f m).', ...
            vehicle, inside{k, 1}, inside{k, 2}, VC.height);
    end
end
end

function tf = hasPath(s, p)
tf = true;
for f = split(string(p), ".")'
    if ~isstruct(s) || ~isfield(s, f), tf = false; return; end
    s = s.(f);
end
end
