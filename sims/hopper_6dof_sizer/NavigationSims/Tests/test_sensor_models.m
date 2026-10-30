%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test_sensor_models.m
%
% Consistency checks between every Src_sim sensor model and its
% Src_filter update model:
%   1. With noise off, the sim measurement equals the filter prediction
%      (innovation = 0) for baro, gyro, GPS, LIDAR and differential GPS
%   2. LIDAR sim ranges equal an independent ray-to-ground intersection
%   3. Every update model's H matches a finite-difference Jacobian taken
%      through state_injection.m, i.e. H = -d h(inject(x, e)) / de
%
% Runs filter_sim_initializer first. Errors if any check fails.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

testDir = fileparts(mfilename('fullpath'));
startDir = cd(fileparts(testDir));
filter_sim_initializer;
rng(3);

% Noise off for every sensor
gps_params.pos_noise_std = zeros(3,1);
gps_params.vel_noise_std = zeros(3,1);
gps_params.multipath = zeros(3,1);
gps_params.bias = zeros(3,1);
imu_params.accel_noise_std = zeros(3,1);
imu_params.gyro_noise_std = zeros(3,1);
baro_params.noise_std = 0;
lidar_params.noise_std_1 = 0;
lidar_params.noise_std_2 = 0;
diffGps_params.noise_std = zeros(3,1);
clear gps_model imu_model barometer_model lidar_model_base diffGps_model

nTrials = 200;
simErr = struct('baro', 0, 'gyro', 0, 'gps', 0, 'lidar', 0, 'diffGps', 0);
lidarGeomErr = 0;
lidarParams = lidar_frames;
for k = 1:nTrials
    x = random_state();
    t = 0.01*k;

    baroBiasTrue = x(20);
    meas = barometer_model(x, baroBiasTrue, t);
    simErr.baro = max(simErr.baro, abs(baro_update_model(x, meas, struct())));

    gyroBiasTrue = x(17:19);
    [~, meas] = imu_model(x, zeros(3,1), zeros(3,1), gyroBiasTrue, t);
    simErr.gyro = max(simErr.gyro, norm(gyro_update_model(x, meas, struct())));

    meas = gps_model(x, t, gps_frame);
    gpsMeas = struct('pos', meas.pos, 'vel', meas.vel);
    gpsParams.r_b = gps_lever_arm(gps_frame);
    simErr.gps = max(simErr.gps, norm(gps_update_model(x, gpsMeas, gpsParams)));

    meas = lidar_system_model(x, t, lidar_frames);
    for i = 1:4
        lidarParams.lidarIndex = i;
        beam = struct('distance', meas.distance(i), 'covariance', 1);
        simErr.lidar = max(simErr.lidar, abs(lidar_update_model(x, beam, lidarParams)));
        lidarGeomErr = max(lidarGeomErr, abs(meas.distance(i) - ray_to_ground(x, lidar_frames, i)));
    end

    meas = diffGps_model(x, t, gps_frame);
    simErr.diffGps = max(simErr.diffGps, norm(diffGps_update_model(x, meas, gps_frame)));
end
fprintf('1. noise-off max |innovation|: baro %.1e, gyro %.1e, gps %.1e, lidar %.1e, diffGps %.1e\n', ...
    simErr.baro, simErr.gyro, simErr.gps, simErr.lidar, simErr.diffGps);
fprintf('2. LIDAR sim vs independent ray intersection: max %.1e m\n', lidarGeomErr);

% Finite-difference Jacobians through state_injection
jacErr = struct('baro', 0, 'gyro', 0, 'gps', 0, 'lidar', 0, 'diffGps', 0);
gpsParams.r_b = gps_lever_arm(gps_frame) + [0.1; 0.05; 0];
for k = 1:50
    x = random_state();
    jacErr.baro = max(jacErr.baro, fd_error(@baro_update_model, x, struct('height', 0), struct()));
    jacErr.gyro = max(jacErr.gyro, fd_error(@gyro_update_model, x, struct('angvel', zeros(3,1)), struct()));
    jacErr.gps = max(jacErr.gps, fd_error(@gps_update_model, x, struct('pos', zeros(3,1), 'vel', zeros(3,1)), gpsParams));
    for i = 1:4
        lidarParams.lidarIndex = i;
        jacErr.lidar = max(jacErr.lidar, fd_error(@lidar_update_model, x, struct('distance', 0), lidarParams));
    end
    jacErr.diffGps = max(jacErr.diffGps, fd_error(@diffGps_update_model, x, struct('basisVec', zeros(3,1)), gps_frame));
end
fprintf('3. max |H - H_fd|: baro %.1e, gyro %.1e, gps %.1e, lidar %.1e, diffGps %.1e\n', ...
    jacErr.baro, jacErr.gyro, jacErr.gps, jacErr.lidar, jacErr.diffGps);

Simulink.data.dictionary.closeAll('-discard');
cd(startDir);

%% Pass/fail
names = fieldnames(simErr);
checks = cell(0, 2);
for i = 1:numel(names)
    checks(end+1,:) = {sprintf('%s sim matches filter (noise off)', names{i}), simErr.(names{i}) < 1e-9}; %#ok<SAGROW>
end
checks(end+1,:) = {'LIDAR sim matches ray intersection', lidarGeomErr < 1e-9};
for i = 1:numel(names)
    checks(end+1,:) = {sprintf('%s Jacobian matches finite difference', names{i}), jacErr.(names{i}) < 1e-5}; %#ok<SAGROW>
end
fprintf('\n');
for i = 1:size(checks,1)
    if checks{i,2}, tag = 'PASS'; else, tag = 'FAIL'; end
    fprintf('  [%s] %s\n', tag, checks{i,1});
end
assert(all([checks{:,2}]), 'test_sensor_models: one or more checks failed');
fprintf('test_sensor_models: all checks passed\n');

%% Helpers
function x = random_state()
% Near-vertical attitude (LIDAR beams must reach the ground), small heights
    tilt = deg2rad(5) * randn(2,1);
    roll = 2*pi*rand();
    e = zeros(19,1);
    e(10:12) = [tilt; roll];
    x = zeros(20,1);
    x(10:13) = [1; 0; 0; 0];
    x = state_injection(x, e);
    x(1:3) = [randn(2,1); -0.5*rand()];          % at most 0.5 m above the pad
    x(4:9) = randn(6,1);
    x(14:20) = 0.01*randn(7,1);
end

function d = ray_to_ground(x, frames, i)
% Distance along beam i from its tip to the ground plane, which is
% frames.heightCOM below the NED origin (positive down)
    C_IB = quaternion_to_rotation(x(10:13));
    A = frames.(sprintf('A%d', i));
    r = frames.(sprintf('r%d', i));
    tip = x(1:3) + C_IB * r;
    beam = C_IB * A(:,1);
    d = (frames.heightCOM - tip(3)) / beam(3);
end

function err = fd_error(model, x, meas, params)
% Max |H - H_fd| with H_fd = -d h(inject(x, e)) / de by central differences
    delta = 1e-6;
    [~, H] = model(x, meas, params);
    Hnum = zeros(size(H));
    for i = 1:19
        e = zeros(19,1);
        e(i) = delta;
        hPlus = model(state_injection(x, e), meas, params);
        hMinus = model(state_injection(x, -e), meas, params);
        Hnum(:,i) = -(hPlus(:) - hMinus(:)) / (2*delta);
    end
    err = max(abs(H - Hnum), [], 'all');
end
