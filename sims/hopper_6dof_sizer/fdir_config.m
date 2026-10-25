function [FDIR_CONFIG, fault_config, fault_map, check_config] = fdir_config()
% FDIR_CONFIG  Single source of truth for FDIR tuning and the fault catalog.
%
% Pure function so the MATLAB Function blocks can call it directly (the
% values fold to constants at code generation) and the tests can call it
% without loading Simulink. Init_FDIR.m pushes the results to the base
% workspace for the model.
%
%   FDIR_CONFIG   rates, sensor layout, thresholds, fault_map
%   fault_config  102 catalog rows: ID, Severity, Response_Tier
%   fault_map     catalog row each of the 16 monitor checks reports against
%   check_config  16 monitor checks: ID, Fault_ID, N_Threshold, M_Window
%
% Persistence runs over the 16 checks, not the 102 catalog rows, because
% several checks share a row and they do not share a threshold: Barometer
% loss has to latch in one frame while Barometer drift still needs 4 of 10,
% and both report against F-04.
%
%#codegen

%% 1. Executive rates
M_WINDOW_SIZE   = 10;     % default rolling window size (M frames)
N_FAIL_THRESH   = 4;      % default failure threshold (N out of M)
SAMPLE_RATE_SEC = 0.01;   % execution rate: 10 ms (100 Hz)

FDIR_CONFIG.M_WINDOW_SIZE   = M_WINDOW_SIZE;
FDIR_CONFIG.N_FAIL_THRESH   = N_FAIL_THRESH;
FDIR_CONFIG.SAMPLE_RATE_SEC = SAMPLE_RATE_SEC;

%% 2. Sensor vector layout
% Sensor_data and fdir_ekf_ref are both this 11-wide stacked vector.
N_SENSORS = 4;
N_CHECKS  = 4;            % loss, drift, saturation, stuck
FDIR_CONFIG.N_SENSORS = N_SENSORS;
FDIR_CONFIG.N_CHECKS  = N_CHECKS;
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

%% 4. Monitor check -> catalog fault row (FaultConditions.xlsx, "Fault List")
% Order is [loss, drift, saturation, stuck] per sensor, sensors in order
% IMU, Barometer, Lidar, GPS. Only F-01..F-08 are sensor rows, so several
% checks share one:
%
%   F-01 IMU loss (no data / stale data)      <- IMU loss, IMU stuck
%   F-02 IMU drift or bias beyond limits      <- IMU drift, IMU saturation
%   F-04 Barometer loss or plume/pressure err <- all four barometer checks
%   F-05 GPS/GNSS dropout                     <- GPS loss, GPS stuck
%   F-06 GNSS multipath / degraded solution   <- GPS drift, GPS saturation
%   F-07 LiDAR altimeter dropout              <- Lidar loss, saturation, stuck
%   F-08 LiDAR disagreement                   <- Lidar drift
%
% Saturation has no dedicated row anywhere: present-but-invalid data is
% folded onto the sensor's degraded-data row. Stale data counts as a
% dropout, which is how the sheet words F-01, F-05 and F-07.
fault_map = uint8([ ...
    1, 2, 2, 1, ...   % IMU       loss, drift, saturation, stuck
    4, 4, 4, 4, ...   % Barometer loss, drift, saturation, stuck
    7, 8, 7, 7, ...   % Lidar     loss, drift, saturation, stuck
    5, 6, 6, 5 ...    % GPS       loss, drift, saturation, stuck
    ]');
FDIR_CONFIG.fault_map = fault_map;

%% 5. Full system lookup table (102 fault conditions)
% Severity and tier come from fdir_fault_catalog_data, which is generated
% from fdir_fault_catalog.csv ("Fault List" sheet of FaultConditions.xlsx)
% by make_fdir_catalog. Do not hand-type these: the previous array was
% typed in the sheet's text-sort order (F-10, F-100, F-101, F-102, F-11,
% ...), which silently shifted every row from F-11 on.
%
% Tier is severity - 1: 1 -> R0, 2 -> R1, 3 -> R2, 4 -> R3, 5 -> R4.
% A severity of 0 is the sheet's N/A (not on the vehicle) and pins to R0.
[severities, tiers] = fdir_fault_catalog_data();

num_faults = numel(severities);
num_checks = N_SENSORS * N_CHECKS;
% Set every FDIR_CONFIG field before the loops below read any of them: code
% generation rejects a field added to a struct that has already been read.
FDIR_CONFIG.NUM_FAULTS = num_faults;
FDIR_CONFIG.NUM_CHECKS = num_checks;

fault_config = repmat(struct( ...
    'ID',            uint8(0), ...
    'Severity',      uint8(0), ...
    'Response_Tier', uint8(0)), num_faults, 1);

for i = 1:num_faults
    fault_config(i).ID            = uint8(i);
    fault_config(i).Severity      = uint8(severities(i));
    fault_config(i).Response_Tier = uint8(tiers(i));
end

%% 6. Per-check persistence thresholds
% Loss and stuck act inside a single frame: a NaN or a frozen reading is
% already in the controller before 4-of-10 persistence (40 ms) could latch,
% and the stuck monitor's own stuck_frames counter is the persistence.
% Drift and saturation keep the default N-of-M.
LOSS = 1; STUCK = 4;      % positions within each sensor's 4 checks

check_config = repmat(struct( ...
    'ID',          uint8(0), ...
    'Fault_ID',    uint8(0), ...
    'N_Threshold', uint8(0), ...
    'M_Window',    uint8(0)), num_checks, 1);

for k = 1:num_checks
    within = mod(k - 1, N_CHECKS) + 1;
    check_config(k).ID       = uint8(k);
    check_config(k).Fault_ID = fault_map(k);
    if within == LOSS || within == STUCK
        check_config(k).N_Threshold = uint8(1);
    else
        check_config(k).N_Threshold = uint8(N_FAIL_THRESH);
    end
    check_config(k).M_Window = uint8(M_WINDOW_SIZE);
end
end
