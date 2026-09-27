function SENS = sensor_params(seed)
%SENSOR_PARAMS  Parameters for the Environment/Sensors measurement models.
%   SENS = sensor_params() returns the struct the sensor blocks read.
%   SENS = sensor_params(seed) sets the noise seed (each sensor derives its
%   own stream from it), e.g. one seed per Monte Carlo run.
%
%   Models follow the Notion pages under Dynamics & GNC > Navigation >
%   Measurement Models. Values marked PLACEHOLDER are datasheet figures or
%   guesses to replace with test data (Allan variance, static logs, CAD).
%
%   Frames: NED with origin on the pad at ground level (ground is D = 0);
%   body x points up the vehicle axis. Lever arms are from the CG, body axes.

if nargin < 1, seed = 1; end
SENS.seed = seed;

%% IMU (ADIS16500 on the flight computer)
% f = (I+Ma)(I+Sa) Rbs f_true + ba + na
% w = (I+Mg)(I+Sg) Rbs w_true + Gg Rbs f_true + bg + ng
% Output: accel m/s^2, gyro deg/s (luna fc_sensors::Imu units).
SENS.imu.Ts        = 0.001;           % 1 kHz
SENS.imu.R_bs      = eye(3);          % body -> sensor axes
SENS.imu.r_b       = [0; 0; 0];       % PLACEHOLDER lever arm from CG (m)
% White noise density -> per-sample sd = density / sqrt(Ts)
SENS.imu.acc_nd    = 0.06 / 60;       % PLACEHOLDER VRW 0.06 m/s/sqrt(hr) -> m/s^2/sqrt(Hz)
SENS.imu.gyro_nd   = deg2rad(0.29) / 60; % PLACEHOLDER ARW 0.29 deg/sqrt(hr) -> rad/s/sqrt(Hz)
% Turn-on bias (drawn once) and in-run Gauss-Markov drift
SENS.imu.acc_b0_sd  = 0.01;           % PLACEHOLDER m/s^2
SENS.imu.gyro_b0_sd = deg2rad(0.05);  % PLACEHOLDER rad/s
SENS.imu.acc_b_sd   = 18e-6 * 9.80665;% PLACEHOLDER bias instability 18 ug
SENS.imu.gyro_b_sd  = deg2rad(8.1 / 3600); % PLACEHOLDER bias instability 8.1 deg/hr
SENS.imu.acc_b_tau  = 100;            % PLACEHOLDER s (from Allan variance)
SENS.imu.gyro_b_tau = 100;            % PLACEHOLDER s
% Scale factor and misalignment (drawn once per run)
SENS.imu.sf_sd      = 1e-3;           % PLACEHOLDER 0.1 %
SENS.imu.misalign_sd = 1e-3;          % PLACEHOLDER rad
SENS.imu.g_sens     = deg2rad(0.01) / 9.80665; % PLACEHOLDER 0.01 deg/s/g, rad/s per m/s^2

%% Magnetometer (LIS2MDL)
% z = A Cbn Bn + beta + eta, Gauss (luna Magnetometer units)
SENS.mag.Ts       = 0.01;             % 100 Hz
SENS.mag.B_n      = [0.224; -0.020; 0.420]; % PLACEHOLDER Earth field at pad, NED, Gauss (NOAA calculator)
SENS.mag.A_sd     = 0.01;             % PLACEHOLDER soft-iron/scale/misalignment spread
SENS.mag.beta_sd  = 0.02;             % PLACEHOLDER hard-iron offset, Gauss
SENS.mag.noise_sd = 0.003;            % PLACEHOLDER 3 mGauss RMS

%% Barometer (MS5611)
% p = p_atm(h) + dp_ground_effect + dp_dynamic + beta + eta, quantized
% Output: pressure Pa, temperature degC (luna Barometer units).
SENS.baro.Ts       = 0.01;            % 100 Hz
SENS.baro.p_pad    = 101325;          % PLACEHOLDER pad pressure, Pa
SENS.baro.T_pad    = 288.15;          % PLACEHOLDER pad temperature, K
SENS.baro.L        = 0.0065;          % lapse rate, K/m
SENS.baro.R        = 287.05;          % J/(kg K)
SENS.baro.g0       = 9.80665;
SENS.baro.r_b      = [0; 0; 0];       % PLACEHOLDER sensor position from CG (m)
SENS.baro.ge_A     = 10;              % PLACEHOLDER ground-effect strength, Pa/m^2
SENS.baro.ge_D     = 3;               % PLACEHOLDER ground-effect height, m
SENS.baro.ge_sign  = 1;               % +1 raises pressure
SENS.baro.k_port   = 0.2;             % PLACEHOLDER fraction of dynamic pressure at the port
SENS.baro.wind_n   = [0; 0; 0];       % wind NED used for airspeed, m/s (set to match the run)
SENS.baro.b0_sd    = 20;              % PLACEHOLDER initial bias, Pa
SENS.baro.b_sd     = 10;              % PLACEHOLDER bias drift size, Pa
SENS.baro.b_tau    = 600;             % PLACEHOLDER s
SENS.baro.noise_sd = 1.2;             % MS5611 OSR 4096 RMS noise ~0.012 mbar
SENS.baro.dq       = 1;               % output resolution, Pa (0.01 mbar)

%% GPS (u-blox ZED-F9P), per the Notion parameter table
% p = p_cg + Cbn' r_b + b_slow + b_fast + eta_p, v = v + Cbn'(w x r_b) + eta_v
% Output: [lat deg, lon deg, alt m, vN, vE, vD m/s, has_fix, num_sats]
SENS.gps.tick     = 0.01;             % block runs at 100 Hz to place fixes and latency
SENS.gps.rate     = 5;                % Hz
SENS.gps.latency  = 0.1;              % PLACEHOLDER s
SENS.gps.mode     = 1;                % 0 standalone, 1 RTK fixed
SENS.gps.r_b      = [1.5; 0; 0];      % PLACEHOLDER antenna lever arm (m), body x up
% Columns: [standalone, RTK]
SENS.gps.sH_slow  = [1.10, 0.006];  SENS.gps.sV_slow  = [2.60, 0.010];
SENS.gps.sH_fast  = [0.60, 0.007];  SENS.gps.sV_fast  = [1.40, 0.011];
SENS.gps.sH_white = [0.10, 0.004];  SENS.gps.sV_white = [0.20, 0.006];
SENS.gps.tau_slow = [1800, 900];      % s (RTK value PLACEHOLDER)
SENS.gps.tau_fast = [60, 60];         % s
SENS.gps.svH      = 0.042;            % velocity noise, m/s
SENS.gps.svD      = 0.068;
SENS.gps.dop_scale = [1; 1];          % [HDOP/HDOPnom; VDOP/VDOPnom]
SENS.gps.p_loss   = [1e-4, 1e-3];     % PLACEHOLDER per-fix dropout chance [engine off, on]
SENS.gps.T_reacq  = 2;                % mean outage, s
SENS.gps.num_sats = 12;               % PLACEHOLDER reported when fixed
SENS.gps.lat0     = 33.7756;          % PLACEHOLDER pad latitude, deg
SENS.gps.lon0     = -84.3963;         % PLACEHOLDER pad longitude, deg
SENS.gps.h0       = 300;              % PLACEHOLDER pad ellipsoid height, m
SENS.gps.dq_deg   = [1e-7, 1e-9];     % NAV-PVT / NAV-HPPOSLLH resolution, deg
SENS.gps.dq_alt   = [1e-3, 1e-4];     % m
SENS.gps.dq_vel   = 1e-3;             % m/s

%% LiDAR (4x Benewake TF-Luna), simplified until the ray-traced model exists
% Range along each beam to the flat ground (D = 0), noise, 1 cm steps,
% random dropouts. Invalid readings are reported as -1.
SENS.lidar.Ts      = 0.01;            % 100 Hz
bot = -1.14;                          % PLACEHOLDER vehicle bottom relative to CG, m (cg_init)
arm = 0.30;                           % PLACEHOLDER radial offset of each sensor, m
SENS.lidar.r_b     = [bot bot bot bot; arm -arm 0 0; 0 0 arm -arm]; % one column per sensor
SENS.lidar.d_b     = repmat([-1; 0; 0], 1, 4);                       % beams point down the body axis
SENS.lidar.r_min   = 0.2;             % m
SENS.lidar.r_max   = 8;               % m
SENS.lidar.sd_near = 0.03;            % m, below 3 m (datasheet +-6 cm taken as 2 sigma)
SENS.lidar.sd_frac = 0.01;            % fraction of range, above 3 m (+-2 % as 2 sigma)
SENS.lidar.dq      = 0.01;            % m
SENS.lidar.p_drop  = 0.01;            % PLACEHOLDER per-sample dropout (plume, dust)
end
