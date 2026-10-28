%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% gmekf_correct_hopper.m
%
% Function to apply sensor measurements to correct the state propagation.
% Hopper version of gmekf_correct.m: differential GPS instead of the
% magnetometer, and the GPS lever arm taken from the antenna frame.
%
% KEEP IN SYNC WITH gmekf_correct.m: pull regularly and copy any new
% sensor handling from gmekf_correct.m into this file.
%
% Based on gmekf_correct.m (David Meisinger, 10/18/26)
%
% Inputs:
%   xnp1: np1 propagated state vector (20x1)
%   Pnp1: np1 propagated covariance (19x19)
%   meas_posvel: GPS measurement struct
%   meas_angvel: gyro measurement struct
%   meas_baro: barometer measurement struct
%   meas_LID: LIDAR measurement struct
%   lidar_frames: LIDAR frame measurement struct
%   meas_diffGps: differential GPS measurement struct
%   gps_frame: GPS antenna positions from COM (gpsFrameBus)
%
% Outputs:
%   xnp1_c: n+1 corrected state (20x1)
%   Pnp1_c: n+1 corrected covariance (19x19)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [xnp1_c, Pnp1_c] = gmekf_correct_hopper( ...
    xnp1, Pnp1, meas_posvel, ...
    meas_angvel, meas_baro, meas_LID, lidar_frames, meas_diffGps, gps_frame)

    % Declare persistent timestamps for last sensor read
    persistent gps_timestamp;
    persistent gyro_timestamp;
    persistent baro_timestamp;
    persistent lid_timestamp;
    persistent diffGps_timestamp;

    % Init persistents if needed
    if isempty(gps_timestamp)
        gps_timestamp = -1;
        gyro_timestamp = -1;
        baro_timestamp = -1;
        lid_timestamp = -1;
        diffGps_timestamp = -1;
    end

    % initialize corrected outputs with propagated values
    xnp1_c = xnp1;
    Pnp1_c = Pnp1;

    % Innovation gate thresholds
    baroGateThreshold = 6.635;      % 99%, 1 DOF
    gyroGateThreshold = 11.345;     % 99%, 3 DOF
    gpsGateThreshold  = 16.812;     % 99%, 6 DOF
    lidarGateThreshold = 6.635;     % 99%, 1 DOF (per beam)
    diffGpsGateThreshold = 11.345;  % 99%, 3 DOF

    % Point the GPS receiver reports, from the antenna frame
    gpsLeverArmBody = gps_lever_arm(gps_frame);

    % Check if each sensor should be incorporated
    gps_incorp = meas_posvel.validity && meas_posvel.timestamp > gps_timestamp;
    gyro_incorp = meas_angvel.validity && meas_angvel.timestamp > gyro_timestamp;
    baro_incorp = meas_baro.validity && meas_baro.timestamp > baro_timestamp;
    lid1_incorp = meas_LID.validity(1) && meas_LID.timestamp > lid_timestamp;
    lid2_incorp = meas_LID.validity(2) && meas_LID.timestamp > lid_timestamp;
    lid3_incorp = meas_LID.validity(3) && meas_LID.timestamp > lid_timestamp;
    lid4_incorp = meas_LID.validity(4) && meas_LID.timestamp > lid_timestamp;
    diffGps_incorp = meas_diffGps.validity && meas_diffGps.timestamp > diffGps_timestamp;

    % Barometer correction
    if baro_incorp

        baroParams.gateThreshold = baroGateThreshold;

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

        gyroParams.gateThreshold = gyroGateThreshold;

        [xnp1_c, Pnp1_c, ~] = update( ...
            xnp1_c, ...
            Pnp1_c, ...
            @gyro_update_model, ...
            meas_angvel, ...
            gyroParams);

        gyro_timestamp = meas_angvel.timestamp;

    end

    if gps_incorp

        gpsParams.gateThreshold = gpsGateThreshold;

        % Antenna lever arm in body coordinates
        gpsParams.r_b = gpsLeverArmBody;

        % TEMPORARY until the GPS bus refactor: gpsBus carries separate
        % position/velocity covariances, update() needs one covariance
        meas_gps = struct('pos', meas_posvel.pos, 'vel', meas_posvel.vel, ...
            'covariance', blkdiag(meas_posvel.pcovariance, meas_posvel.vcovariance));

        [xnp1_c, Pnp1_c, ~] = update( ...
            xnp1_c, ...
            Pnp1_c, ...
            @gps_update_model, ...
            meas_gps, ...
            gpsParams);

        gps_timestamp = meas_posvel.timestamp;

    end

    % LIDAR correction
    lidarParams = lidar_frames;
    lidarParams.gateThreshold = lidarGateThreshold;
    lidarParams.lidarIndex = 0;
    if lid1_incorp

        meas_lid1.distance = meas_LID.distance(1);
        meas_lid1.covariance = meas_LID.covariance(1);
        lidarParams.lidarIndex = 1;

        [xnp1_c, Pnp1_c, ~] = update( ...
            xnp1_c, ...
            Pnp1_c, ...
            @lidar_update_model, ...
            meas_lid1, ...
            lidarParams);

        lid_timestamp = meas_LID.timestamp;

    end
    if lid2_incorp

        meas_lid2.distance = meas_LID.distance(2);
        meas_lid2.covariance = meas_LID.covariance(2);
        lidarParams.lidarIndex = 2;

        [xnp1_c, Pnp1_c, ~] = update( ...
            xnp1_c, ...
            Pnp1_c, ...
            @lidar_update_model, ...
            meas_lid2, ...
            lidarParams);

        lid_timestamp = meas_LID.timestamp;

    end
    if lid3_incorp

        meas_lid3.distance = meas_LID.distance(3);
        meas_lid3.covariance = meas_LID.covariance(3);
        lidarParams.lidarIndex = 3;

        [xnp1_c, Pnp1_c, ~] = update( ...
            xnp1_c, ...
            Pnp1_c, ...
            @lidar_update_model, ...
            meas_lid3, ...
            lidarParams);

        lid_timestamp = meas_LID.timestamp;

    end
    if lid4_incorp

        meas_lid4.distance = meas_LID.distance(4);
        meas_lid4.covariance = meas_LID.covariance(4);
        lidarParams.lidarIndex = 4;

        [xnp1_c, Pnp1_c, ~] = update( ...
            xnp1_c, ...
            Pnp1_c, ...
            @lidar_update_model, ...
            meas_lid4, ...
            lidarParams);

        lid_timestamp = meas_LID.timestamp;

    end

    % Differential GPS correction (antenna 1 -> antenna 2 baseline in NED)
    if diffGps_incorp

        diffGpsParams.gateThreshold = diffGpsGateThreshold;

        % Antenna positions in body coordinates
        diffGpsParams.r1 = gps_frame.r1;
        diffGpsParams.r2 = gps_frame.r2;

        [xnp1_c, Pnp1_c, ~] = update( ...
            xnp1_c, ...
            Pnp1_c, ...
            @diffGps_update_model, ...
            meas_diffGps, ...
            diffGpsParams);

        diffGps_timestamp = meas_diffGps.timestamp;

    end

end
