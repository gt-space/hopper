%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% gmekf_correct_hopper.m
%
% Function to apply sensor measurements to correct the state propagation.
% Hopper version of gmekf_correct.m: differential GPS instead of the
% magnetometer, and the GPS lever arm taken from the antenna frame.
%
% Based on gmekf_correct.m (David Meisinger, 10/18/26)
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
%   meas_diffGps: differential GPS measurement struct
%   gps_frame: GPS antenna positions from COM (gpsFrameBus)
%
% Outputs:
%   xnp1_c: n+1 corrected state (20x1)
%   Pnp1_c: n+1 corrected covariance (19x19)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [xnp1_c, Pnp1_c] = gmekf_correct_hopper( ...
    xnp1, Pnp1, aux, meas_posvel, ...
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

    % LIDAR is not implemented yet (lidar_frames unused until then).
    % Keep the 4x1 validity size so the bus type doesn't change.
    meas_LID.validity(:) = false;
    meas_LID.timestamp = -1;

    % initialize corrected outputs with propagated values
    xnp1_c = xnp1;
    Pnp1_c = Pnp1;

    aux.baroGateThreshold = 6.635;      % 99%, 1 DOF
    aux.gyroGateThreshold = 11.345;     % 99%, 3 DOF
    aux.gpsGateThreshold  = 16.812;     % 99%, 6 DOF
    aux.diffGpsGateThreshold = 11.345;  % 99%, 3 DOF

    % GPS position is the antenna midpoint, matching gps_model.m
    aux.gpsLeverArmBody = (gps_frame.r1(:) + gps_frame.r2(:))/2;

    %{
        Measurement adapters: the *_update_model functions and update.m
        expect different field names than the dictionary buses (update.m
        reads R, the buses carry covariance; gyro_update_model reads
        angularVelocity, the bus carries angvel). Build the structs they
        expect here. Remove once the names are reconciled.
    %}
    baroMeas = struct('height', meas_baro.height, 'R', meas_baro.covariance);
    gyroMeas = struct('angularVelocity', meas_angvel.angvel, 'R', meas_angvel.covariance);
    gpsMeas = struct('pos', meas_posvel.pos, 'vel', meas_posvel.vel, ...
        'R', blkdiag(meas_posvel.pcovariance, meas_posvel.vcovariance));
    diffGpsMeas = struct('basisVec', meas_diffGps.basisVec, 'R', meas_diffGps.covariance);

    % Check if each sensor should be incorporated
    gps_incorp = meas_posvel.validity && meas_posvel.timestamp > gps_timestamp;
    gyro_incorp = meas_angvel.validity && meas_angvel.timestamp > gyro_timestamp;
    baro_incorp = meas_baro.validity && meas_baro.timestamp > baro_timestamp;
    lid_incorp = all(meas_LID.validity) && meas_LID.timestamp > lid_timestamp;
    diffGps_incorp = meas_diffGps.validity && meas_diffGps.timestamp > diffGps_timestamp;


    % Barometer correction
    if baro_incorp

        baroParams.gateThreshold = aux.baroGateThreshold;

        % Perform complete GMEKF update
        [xnp1_c, Pnp1_c, ~] = update( ...
            xnp1_c, ...
            Pnp1_c, ...
            @baro_update_model, ...
            baroMeas, ...
            baroParams);

        baro_timestamp = meas_baro.timestamp; % consume measurement even if gated out
    end

    if gyro_incorp

        gyroParams.gateThreshold = aux.gyroGateThreshold;

        [xnp1_c, Pnp1_c, ~] = update( ...
            xnp1_c, ...
            Pnp1_c, ...
            @gyro_update_model, ...
            gyroMeas, ...
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
            gpsMeas, ...
            gpsParams);

        gps_timestamp = meas_posvel.timestamp;

    end

    % Differential GPS correction (antenna baseline in NED)
    if diffGps_incorp

        diffGpsParams.gateThreshold = aux.diffGpsGateThreshold;

        % Antenna positions in body coordinates
        diffGpsParams.r1 = gps_frame.r1;
        diffGpsParams.r2 = gps_frame.r2;

        [xnp1_c, Pnp1_c, ~] = update( ...
            xnp1_c, ...
            Pnp1_c, ...
            @diffGps_update_model, ...
            diffGpsMeas, ...
            diffGpsParams);

        diffGps_timestamp = meas_diffGps.timestamp;

    end

    % TODO: LIDAR correction

end
