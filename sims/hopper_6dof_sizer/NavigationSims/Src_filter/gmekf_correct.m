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
function [xnp1_c, Pnp1_c] = gmekf_correct(xnp1, Pnp1, aux, meas_posvel, ...
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

    % Check if each sensor should be incorporated
    gps_incorp = meas_posvel.validity && meas_posvel.timestamp ~= gps_timestamp;
    gyro_incorp = meas_angvel.validity && meas_angvel.timestamp ~= gyro_timestamp;
    baro_incorp = meas_baro.validity && meas_barovel.timestamp ~= baro_timestamp;
    lid_incorp = meas_LID.validity && meas_LID.timestamp ~= lid_timestamp;
    mag_incorp = meas_mag.validity && meas_mag.timestamp ~= mag_timestamp;

    % Update the state estimate

end