function S = sim_sensor_params(IN)
% SIM_SENSOR_PARAMS  Sensor noise, rates and geometry for the Sensors block.
%
%   S = sim_sensor_params(IN)   IN.sensors overrides any field
%
% Noise values follow navigation's filter_sim_initializer (IMU, baro, GPS);
% the lidar noise is FDIR's assumed 0.02 m (navigation's lidar model has
% none yet) and must not be zero, or a vehicle sitting on the pad trips
% the stuck checks. Rates follow the GNC PDR (slide 87) where it gives
% them, capped at the 100 Hz FDIR frame.
%
% Mounting geometry (lidar azimuth, cant, radius, mount_x) is vehicle data:
% it comes from the vehicle config's sensors section through IN.sensors
% (vehicles/<vehicle>/<vehicle>_config.m), in the sim's body axes.

S.frame = 0.01;                                % s, FDIR / sensor frame (100 Hz)

S.imu.period      = 0.01;                      % s (IMU is 1 kHz on the vehicle)
S.imu.accel_sigma = [0.05; 0.05; 0.05];        % m/s^2
S.imu.gyro_sigma  = [0.005; 0.005; 0.005];     % rad/s
S.imu.accel_fs    = 16 * 9.80665;              % m/s^2, +/-16 g
S.imu.gyro_fs     = deg2rad(2000);             % rad/s

S.baro.period = 0.05;                          % s
S.baro.sigma  = 0.1;                           % m

S.lidar.period   = 0.01;                       % s
S.lidar.sigma    = 0.02;                       % m
S.lidar.r_min    = 0.05;                       % m, blind zone (masked on the pad)
S.lidar.r_max    = 40;                         % m

S.gps.period = 0.2;                            % s (GNSS 5 Hz)
S.gps.sigma  = [0.02; 0.02; 0.05];             % m
S.gps.fix    = 3;                              % 3D
S.gps.sats   = 12;

S.mag.period = 0.02;                           % s
S.mag.sigma  = 0.5;                            % uT
S.mag.field_ned = [22.6; -2.5; 42.0];          % uT, Earth field near Atlanta (WMM, approx.)

S.seed = 42;                                   % noise seed; Monte Carlo can vary it

if isfield(IN, 'sensors')
    S = mergeStruct(S, IN.sensors);
end
if ~all(isfield(S.lidar, {'azimuth', 'cant', 'radius', 'mount_x'}))
    error('sim_sensor_params:geometry', 'Lidar mounting comes from the vehicle config (VC.sensors.lidar); IN.sensors has none.');
end
end

function a = mergeStruct(a, b)
for f = fieldnames(b)'
    if isstruct(b.(f{1})) && isfield(a, f{1}) && isstruct(a.(f{1}))
        a.(f{1}) = mergeStruct(a.(f{1}), b.(f{1}));
    else
        a.(f{1}) = b.(f{1});
    end
end
end
