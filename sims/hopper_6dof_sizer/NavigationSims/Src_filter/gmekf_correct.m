%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% gmekf_correct.m
% 
% Function to apply sensor measurements to correct the state propagation.
% 
% Author: David Meisinger
% Date: 10/18/26
% 
% Inputs:
%   xnp1: np1 propagated state vector (20x1)
%   Pnp1: np1 propagated covariance (19x19)
%   aux: Struct for passing common Navigation params
%   meas_posvel: GPS measurement struct
%   meas_angvel: gyro measurement struct
%   meas_baro: barometer measurement struct
%   meas_LID: LIDAR measurement struct
%   lidar_frames: LIDAR frame measurement struct
%   meas_mag: magnetometer measurement struct
% 
% Outputs:
%   xnp1_c: n+1 corrected state (20x1)
%   Pnp1_c: n+1 corrected covariance (19x19)
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [xnp1_c, Pnp1_c] = gmekf_correct( ...
    xnp1, Pnp1, aux, meas_posvel, ...
    meas_angvel, meas_baro, meas_LID, lidar_frames, meas_mag)

    % Declare persistent timestamps for last sensor read
    persistent gps_timestamp;
    persistent gyro_timestamp;
    persistent baro_timestamp;
    persistent lid_timestamp;
    persistent mag_timestamp;

    % Init persistents if needed
    if isempty(gps_timestamp)
        gps_timestamp = -1;
        gyro_timestamp = -1;
        baro_timestamp = -1;
        lid_timestamp = -1;
        mag_timestamp = -1;
    end

    % LIDAR is not implemented yet
    meas_LID.validity = false;
    meas_LID.timestamp = -1;

    % Magnetometer is not implemented yet
    meas_mag.validity = false;
    meas_mag.timestamp = -1;

    % Unused LIDAR frame data
    lidar_frames = struct();

    % initialize corrected outputs with propagated values
    xnp1_c = xnp1;
    Pnp1_c = Pnp1;

    aux.baroGateThreshold = 6.635;   % 99%, 1 DOF
    aux.gyroGateThreshold = 11.345;  % 99%, 3 DOF
    aux.gpsGateThreshold  = 16.812;  % 99%, 6 DOF

    aux.gpsLeverArmBody = [0.1; 0; -0.2];

    % Check if each sensor should be incorporated
    gps_incorp = meas_posvel.validity && meas_posvel.timestamp > gps_timestamp;
    gyro_incorp = meas_angvel.validity && meas_angvel.timestamp > gyro_timestamp;
    baro_incorp = meas_baro.validity && meas_baro.timestamp > baro_timestamp;
    lid_incorp = any(meas_LID.validity) && meas_LID.timestamp > lid_timestamp;
    mag_incorp = meas_mag.validity && meas_mag.timestamp > mag_timestamp;


    % Barometer correction
    if baro_incorp

        baroParams.gateThreshold = aux.baroGateThreshold;

        % Perform complete GMEKF update
        [xnp1_c, Pnp1_c, ~] = update( ...
            xnp1_c, ...
            Pnp1_c, ...
            @baro_update_model, ...
            meas_baro, ...
            baroParams);

        baro_timestamp = meas_baro.timestamp; % consume measurement even if gated out
    end

    if gyro_incorp

        gyroParams.gateThreshold = aux.gyroGateThreshold;

        [xnp1_c, Pnp1_c, ~] = update( ...
            xnp1_c, ...
            Pnp1_c, ...
            @gyro_update_model, ...
            meas_angvel, ...
            gyroParams);

        gyro_timestamp = meas_angvel.timestamp;

    end

    if gps_incorp

        gpsParams.gateThreshold = aux.gpsGateThreshold;

        % Antenna lever arm in body coordinates
        gpsParams.r_b = aux.gpsLeverArmBody;

        [xnp1_c, Pnp1_c, ~] = update( ...
            xnp1_c, ...
            Pnp1_c, ...
            @gps_update_model, ...
            meas_posvel, ...
            gpsParams);

        gps_timestamp = meas_posvel.timestamp;

    end

    % TODO: LIDAR correction
    % TODO: Magnetometer correction
    % TODO: Differential GPS

end