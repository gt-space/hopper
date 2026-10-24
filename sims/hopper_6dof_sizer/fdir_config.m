function [FDIR_CONFIG, fault_config, fault_map] = fdir_config()
% FDIR_CONFIG  Single source of truth for FDIR tuning and the fault catalog.
%
% Pure function so the MATLAB Function blocks can call it directly (the
% values fold to constants at code generation) and the tests can call it
% without loading Simulink. Init_FDIR.m pushes the results to the base
% workspace for the model.
%
% [FDIR_CONFIG, fault_config, fault_map] = fdir_config()
%
%#codegen

%% 1. Executive rates
FDIR_CONFIG.M_WINDOW_SIZE  = 10;      % default rolling window size (M frames)
FDIR_CONFIG.N_FAIL_THRESH  = 4;       % default failure threshold (N out of M)
FDIR_CONFIG.SAMPLE_RATE_SEC = 0.01;   % execution rate: 10 ms (100 Hz)

%% 2. Sensor vector layout
% Sensor_data and fdir_ekf_ref are both this 11-wide stacked vector.
FDIR_CONFIG.N_SENSORS = 4;
FDIR_CONFIG.N_CHECKS  = 4;            % loss, drift, saturation, stuck
FDIR_CONFIG.WIDTH     = 11;

%% 3. Per-sensor thresholds
% drift is a per-channel band on |measurement - EKF reference|, sized at
% ~5 sigma of that channel's noise. No sensor noise model lives in this repo
% yet, so the sigmas below are the assumed values and are called out per
% channel; revisit once the EKF noise model lands.

% IMU: 3 accel (m/s^2) then 3 gyro (rad/s).
% sigma assumed 0.02 m/s^2 per accel axis and 0.002 rad/s per gyro axis
% (industrial MEMS IMU, 100 Hz, in-run bias excluded).
s.IMU.idx          = (1:6);
s.IMU.drift        = [0.10; 0.10; 0.10; 0.010; 0.010; 0.010];
s.IMU.range_lo     = [-160; -160; -160; -34.9; -34.9; -34.9];   % +/-16 g, +/-2000 deg/s
s.IMU.range_hi     = [ 160;  160;  160;  34.9;  34.9;  34.9];
s.IMU.stuck_frames = 20;                                        % 0.2 s

% Barometer: pressure altitude (m).
% sigma assumed 0.3 m (MS5611-class, oversampled, 100 Hz).
s.Barometer.idx          = 7;
s.Barometer.drift        = 1.5;
s.Barometer.range_lo     = -500;
s.Barometer.range_hi     = 9000;
s.Barometer.stuck_frames = 20;

% Lidar: ground range (m).
% sigma assumed 0.02 m (single-point ToF rangefinder at close range).
s.Lidar.idx          = 8;
s.Lidar.drift        = 0.10;
s.Lidar.range_lo     = 0;
s.Lidar.range_hi     = 120;
s.Lidar.stuck_frames = 20;

% GPS: NED position (m).
% sigma assumed 1.5 m horizontal and 3.0 m vertical (SBAS-corrected L1).
s.GPS.idx          = (9:11);
s.GPS.drift        = [7.5; 7.5; 15.0];
s.GPS.range_lo     = [-5e4; -5e4; -2e4];
s.GPS.range_hi     = [ 5e4;  5e4;  5e2];
s.GPS.stuck_frames = 20;

FDIR_CONFIG.sensors = s;

%% 4. Monitor check -> catalog fault ID
% TODO confirm against FaultConditions.xlsx. FaultConditions.xlsx is not in
% the repo, so this block is a placeholder: the 16 monitor checks are mapped
% onto catalog rows F-01..F-16 in sensor order. The stuck checks in
% particular may need new catalog rows, since a frozen-sensor condition may
% not exist in the current 102.
%
% Order is [loss, drift, saturation, stuck] per sensor,
% sensors in order IMU, Barometer, Lidar, GPS.
fault_map = uint8([ ...
     1,  2,  3,  4, ...   % IMU       loss, drift, saturation, stuck
     5,  6,  7,  8, ...   % Barometer loss, drift, saturation, stuck
     9, 10, 11, 12, ...   % Lidar     loss, drift, saturation, stuck
    13, 14, 15, 16 ...    % GPS       loss, drift, saturation, stuck
    ]');
FDIR_CONFIG.fault_map = fault_map;

%% 5. Full system lookup table (102 fault conditions)
% Direct mapping from FaultConditions.xlsx:
% Severity 1 -> R0 (Nominal/Mask)
% Severity 2 -> R1 (Degrade / Mask)
% Severity 3 -> R2 (Abort to Land)
% Severity 4 -> R3 (Nav-Degraded Abort)
% Severity 5 -> R4 (Cutoff & Safing)
severities = [ ...
    3, 3, 2, 3, 3, 2, 3, 3, 4, 3, 5, 5, 5, 3, 1, 5, 4, 5, 4, 4, ... % F-01..F-20
    4, 3, 4, 4, 4, 5, 3, 2, 2, 2, 2, 5, 2, 2, 4, 2, 4, 3, 5, 5, ... % F-21..F-40
    5, 5, 3, 5, 5, 4, 3, 4, 5, 2, 2, 2, 1, 1, 2, 2, 2, 2, 3, 3, ... % F-41..F-60
    3, 5, 2, 3, 3, 3, 3, 2, 2, 2, 5, 2, 5, 5, 1, 5, 2, 4, 2, 2, ... % F-61..F-80
    3, 2, 3, 1, 5, 3, 1, 1, 2, 2, 3, 4, 5, 5, 3, 4, 4, 3, 3, 5, ... % F-81..F-100
    5, 5 ...                                                        % F-101..F-102
    ];

num_faults = numel(severities);
% Set every FDIR_CONFIG field before the loop below reads any of them: code
% generation rejects a field added to a struct that has already been read.
FDIR_CONFIG.NUM_FAULTS = num_faults;

% Faults that have to act inside a single frame get N = 1: a NaN or an Inf
% propagates into the controller on the next frame, well before 4-of-10
% persistence (40 ms) could latch. The stuck checks are also N = 1 because
% the monitor's own stuck_frames counter already provides the persistence.
fast_path = false(1, num_faults);
fast_path(fault_map(1:4:end)) = true;   % loss checks
fast_path(fault_map(4:4:end)) = true;   % stuck checks

fault_config = repmat(struct( ...
    'ID',            uint8(0), ...
    'Severity',      uint8(0), ...
    'N_Threshold',   uint8(0), ...
    'M_Window',      uint8(0), ...
    'Response_Tier', uint8(0)), num_faults, 1);

for i = 1:num_faults
    fault_config(i).ID       = uint8(i);
    fault_config(i).Severity = uint8(severities(i));

    if fast_path(i)
        fault_config(i).N_Threshold = uint8(1);
    else
        fault_config(i).N_Threshold = uint8(FDIR_CONFIG.N_FAIL_THRESH);
    end
    fault_config(i).M_Window = uint8(FDIR_CONFIG.M_WINDOW_SIZE);

    % Response tier derived from severity (1->R0, 2->R1, 3->R2, 4->R3, 5->R4)
    fault_config(i).Response_Tier = uint8(max(0, severities(i) - 1));
end
end
