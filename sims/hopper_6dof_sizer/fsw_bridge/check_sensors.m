function check_sensors()
%CHECK_SENSORS  Sanity-check the sensor models against the true state.
%   Runs hopper_6dof_NED_v2_fswBridge and compares each measurement with the
%   truth it is built from, printing error statistics next to what the
%   parameters in sensor_params imply. The model is never saved.

here = fileparts(mfilename('fullpath'));
cd(fileparts(here));
if ~evalin('base', 'exist(''IN'', ''var'')'), evalin('base', 'sim_setup_cached'); end
if ~evalin('base', 'exist(''SENS'', ''var'')'), assignin('base', 'SENS', sensor_params()); end
S = evalin('base', 'SENS');

mdl = 'hopper_6dof_NED_v2_fswBridge';
env = [mdl '/Environment'];
load_system(mdl);
% Log the true state and the plant signals the sensors read
sen = [env '/Sensors'];
names = {'F_b', 'W_b', 'mass', 'w_b', 'p_n', 'v_n', 'C_bn', 'thrust'};
for i = 1:numel(names)
    lh = get_param([sen '/' names{i}], 'LineHandles');
    src = get_param(lh.Outport, 'SrcPortHandle');
    set_param(src, 'DataLogging', 'on', 'DataLoggingNameMode', 'Custom', 'DataLoggingName', ['true_' names{i}]);
end
set_param(mdl, 'SignalLogging', 'on', 'SignalLoggingName', 'logsout');
out = sim(Simulink.SimulationInput(mdl));

tr = @(n) out.logsout.get(['true_' n]).Values;
at = @(ts, t) interp1(ts.Time, rows(ts.Data, numel(ts.Time)), t, 'previous');
fprintf('Flight: %.2f s\n', out.get('z').Time(end));

%% IMU (compare where the truth is steady within a sample)
imu = out.get('imu_meas'); t = imu.Time; I = rows(imu.Data, numel(t));
F = at(tr('F_b'), t); W = at(tr('W_b'), t); m = at(tr('mass'), t); w = at(tr('w_b'), t);
f_true = (F - W) ./ m;
ea = I(:, 1:3) - f_true;
eg = I(:, 4:6) - w * 180 / pi;
sa = S.imu.acc_nd / sqrt(S.imu.Ts); sg = S.imu.gyro_nd / sqrt(S.imu.Ts) * 180 / pi;
fprintf('\nIMU (%d samples)\n', numel(t));
fprintf('  accel error  mean %s  std %s m/s^2   (white noise sd %.3g)\n', v3(mean(ea)), v3(std(ea)), sa);
fprintf('  gyro  error  mean %s  std %s deg/s   (white noise sd %.3g)\n', v3(mean(eg)), v3(std(eg)), sg);
fprintf('  accel x at t=5 s: %.3f m/s^2 (thrust/mass %.3f)\n', I(find(t >= 5, 1), 1), ...
    at(tr('thrust'), 5) / at(tr('mass'), 5));

%% Magnetometer
mg = out.get('mag_meas'); M = rows(mg.Data, numel(mg.Time));
fprintf('\nMagnetometer: |B| mean %.4f G (true %.4f), per-axis noise-free spread from A, beta\n', ...
    mean(vecnorm(M, 2, 2)), norm(S.mag.B_n));

%% Barometer: invert the standard atmosphere and compare heights
ba = out.get('baro_meas'); tb = ba.Time; B = rows(ba.Data, numel(tb));
P = S.baro;
h_meas = P.T_pad / P.L * (1 - (B(:, 1) / P.p_pad).^(P.L * P.R / P.g0));
p_n = at(tr('p_n'), tb);
h_true = -p_n(:, 3);
eh = h_meas - h_true;
fprintf('\nBarometer: height error mean %.2f m, std %.2f m, max %.2f m (includes ground effect and bias)\n', ...
    mean(eh), std(eh), max(abs(eh)));

%% GPS: back to NED and compare with the antenna position
gp = out.get('gps_meas'); tg = gp.Time; G = rows(gp.Data, numel(tg));
fixed = G(:, 7) == 1 & [true; any(diff(G(:, 1:3)) ~= 0, 2)];
Q = S.gps; lat0 = Q.lat0 * pi / 180; a = 6378137; e2 = 0.00669438;
N0 = a / sqrt(1 - e2 * sin(lat0)^2); RM = a * (1 - e2) / (1 - e2 * sin(lat0)^2)^1.5;
pN = (G(:, 1) - Q.lat0) * pi / 180 * (RM + Q.h0);
pE = (G(:, 2) - Q.lon0) * pi / 180 * (N0 + Q.h0) * cos(lat0);
pD = Q.h0 - G(:, 3);
% Each fix describes the state `latency` before it is output
tf = tg(fixed) - Q.latency;
C = at(tr('C_bn'), tf);
pt = at(tr('p_n'), tf);
ant = zeros(numel(tf), 3);
for k = 1:numel(tf)
    ant(k, :) = pt(k, :) + (reshape(C(k, :), 3, 3)' * Q.r_b)';
end
ep = [pN(fixed) pE(fixed) pD(fixed)] - ant;
mode = {'standalone', 'RTK'};
fprintf('\nGPS (%s, %d fixes, %d without fix): position error mean %s std %s m\n', ...
    mode{Q.mode + 1}, nnz(fixed), nnz(G(:, 7) == 0), v3(mean(ep)), v3(std(ep)));

%% LiDAR: compare valid ranges with the geometric height of each sensor
li = out.get('lidar_meas'); tl = li.Time; L = rows(li.Data, numel(tl));
p_n = at(tr('p_n'), tl); C = at(tr('C_bn'), tl);
err = []; nvalid = 0;
for k = 1:numel(tl)
    Cnb = reshape(C(k, :), 3, 3)';
    for i = 1:4
        if L(k, i) < 0, continue; end
        pos = p_n(k, :)' + Cnb * S.lidar.r_b(:, i);
        d = Cnb * S.lidar.d_b(:, i);
        err(end + 1) = L(k, i) - (-pos(3)) / d(3); %#ok<AGROW>
        nvalid = nvalid + 1;
    end
end
fprintf('\nLiDAR: %d valid of %d readings, range error mean %.3f m std %.3f m\n', ...
    nvalid, 4 * numel(tl), mean(err), std(err));
end

function A = rows(D, nT)
if ndims(D) == 3
    A = reshape(permute(D, [3 1 2]), nT, []);
elseif size(D, 1) == nT
    A = D;
else
    A = D.';
end
end

function s = v3(x)
s = sprintf('[%s]', strjoin(arrayfun(@(v) sprintf('%.3g', v), x, 'UniformOutput', false), ' '));
end
