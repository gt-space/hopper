function S = sim_sensor_params(IN)
% SIM_SENSOR_PARAMS  Sensor noise, rates and geometry for the Sensors block.
%
%   S = sim_sensor_params(IN)   IN.sensors overrides any field
%   S = sim_sensor_params()    shared defaults, including IMU calibration
%
% IMU calibration values live here once: sample standard deviations and
% continuous-time bias random-walk PSDs. The GMEKF uses these same values
% for its process covariance. Lidar noise is FDIR's assumed 0.02 m
% (navigation's lidar model has none yet) and must not be zero, or a vehicle
% sitting on the pad trips the stuck checks. Rates follow the GNC PDR
% (slide 87) where it gives them, capped at the 100 Hz FDIR frame.
%
% Mounting geometry (lidar azimuth, cant, radius, mount_x; GPS antennas) and
% the IMU rate are vehicle data:
% it comes from the vehicle config's sensors section through IN.sensors
% (vehicles/<vehicle>/<vehicle>_config.m), in the sim's body axes.

if nargin < 1, IN = []; end

S.frame = 0.01;                                % s, FDIR / sensor frame (100 Hz)

S.imu.period      = 1/500;                     % s; vehicle config may override
S.imu.noise_reference_period = 1/500;           % s, calibration sample rate
S.imu.accel_sigma = [0.00759; 0.00843; 0.0119]; % m/s^2, calibrated sample SD at 500 Hz
S.imu.gyro_sigma  = [0.00137; 0.00126; 0.000918]; % rad/s, calibrated sample SD at 500 Hz
S.imu.accel_bias_rw_psd = [9.09e-5; 7.77e-5; 7.16e-5]; % (m/s^2)^2/s, Allan K estimates
S.imu.gyro_bias_rw_psd = [3.52e-6; 2.97e-6; 4.36e-6]; % (rad/s)^2/s, Allan K estimates
S.imu.accel_fs    = 16 * 9.80665;              % m/s^2, +/-16 g
S.imu.gyro_fs     = deg2rad(2000);             % rad/s

S.baro.period = 0.05;                          % s
S.baro.sigma  = 0.1;                           % m
S.baro.bias_noise_psd = 0.001;                 % m^2/s

S.lidar.period   = 0.01;                       % s
S.lidar.sigma    = 0.02;                       % m
S.lidar.r_min    = 0.05;                       % m, blind zone (masked on the pad)
S.lidar.r_max    = 40;                         % m

S.gps.period = 0.2;                            % s (GNSS 5 Hz)
S.gps.sigma  = [1.5; 1.61; 4.42];              % m
S.gps.vel_sigma = [0.029; 0.032; 0.023];       % m/s, receiver placeholder
S.gps.fix    = 3;                              % 3D
S.gps.sats   = 12;

S.dgps.period = 0.1;                           % s, differential GPS (dual antenna) baseline
S.dgps.sigma  = [0.007; 0.007; 0.012];         % m per NED axis (nav team's placeholder)

S.mag.period = 0.02;                           % s
S.mag.sigma  = 0.5;                            % uT
S.mag.field_ned = [22.6; -2.5; 42.0];          % uT, Earth field near Atlanta (WMM, approx.)

S.seed = 42;                                   % noise seed; Monte Carlo can vary it

if isfield(IN, 'sensors')
    S = mergeStruct(S, IN.sensors);
end
if ~isscalar(S.imu.period) || ~isfinite(S.imu.period) || S.imu.period <= 0
    error('sim_sensor_params:imuPeriod', 'IMU sample period must be positive.');
end
sampleScale = sqrt(S.imu.noise_reference_period / S.imu.period);
S.imu.accel_sigma = S.imu.accel_sigma(:) * sampleScale;
S.imu.gyro_sigma = S.imu.gyro_sigma(:) * sampleScale;
if ~isempty(IN) && (~all(isfield(S.lidar, {'azimuth', 'cant', 'radius', 'mount_x'})) || ~isfield(S.gps, 'antenna'))
    error('sim_sensor_params:geometry', 'Lidar and GPS antenna mounting come from the vehicle config (VC.sensors); IN.sensors has none.');
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
