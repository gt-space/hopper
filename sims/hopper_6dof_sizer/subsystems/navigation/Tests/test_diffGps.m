%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test_diffGps.m
%
% Checks for the differential GPS sim model (diffGps_model.m) and filter
% model (diffGps_update_model.m):
%   1. With noise off, the sim measurement equals the filter prediction
%   2. H matches a finite-difference Jacobian taken through
%      state_injection.m, i.e. H = -d h(inject(x, e)) / de
%   3. Sim noise matches the reported covariance
%   4. One update() call removes a roll error about the vertical axis and
%      leaves rotation about the baseline axis unobserved
%
% Runs filter_sim_initializer first. Errors if any check fails.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

testDir = fileparts(mfilename('fullpath'));
startDir = cd(fileparts(testDir));
filter_sim_initializer;
rng(1);

frame = gps_frame;
baselineBody = frame.r2 - frame.r1;
meas0.basisVec = zeros(3,1);    % zero measurement so innovation = h(x)

%% 1. Noise-off sim truth equals filter prediction
diffGps_params_noisy = diffGps_params;
diffGps_params.noise_std = zeros(3,1);
clear diffGps_model             % reset persistent params

truthErr = 0;
lengthErr = 0;
for k = 1:500
    x = random_state();
    meas = diffGps_model(x, 0.1*k, frame);
    innovation = diffGps_update_model(x, meas, frame);
    truthErr = max(truthErr, norm(innovation));
    lengthErr = max(lengthErr, abs(norm(meas.basisVec) - norm(baselineBody)));
end
fprintf('1. noise-off |predicted - measured| max = %.2e m, baseline length err max = %.2e m\n', ...
    truthErr, lengthErr);

%% 2. Finite-difference Jacobian through state_injection
delta = 1e-6;
jacErr = 0;
flippedErr = 0;
for trial = 1:200
    x = random_state();
    [~, H] = diffGps_update_model(x, meas0, frame);
    Hnum = zeros(3,19);
    for i = 1:19
        e = zeros(19,1);
        e(i) = delta;
        hPlus = diffGps_update_model(state_injection(x, e), meas0, frame);
        hMinus = diffGps_update_model(state_injection(x, -e), meas0, frame);
        Hnum(:,i) = -(hPlus - hMinus) / (2*delta);
    end
    jacErr = max(jacErr, max(abs(H - Hnum), [], 'all'));
    Hflip = H;
    Hflip(:,10:12) = -Hflip(:,10:12);
    flippedErr = max(flippedErr, max(abs(Hflip - Hnum), [], 'all'));
end
fprintf('2. max |H - H_fd| = %.2e   (control, attitude sign flipped: %.2e)\n', jacErr, flippedErr);

% Same check on the GPS position rows, for reference only
gpsJacErr = 0;
for trial = 1:50
    x = random_state();
    gpsMeas.pos = zeros(3,1);
    gpsMeas.vel = zeros(3,1);
    gpsParams.r_b = (frame.r1 + frame.r2)/2 + [0.1; 0.05; 0];
    [~, H] = gps_update_model(x, gpsMeas, gpsParams);
    Hnum = zeros(6,19);
    for i = 1:19
        e = zeros(19,1);
        e(i) = delta;
        hPlus = gps_update_model(state_injection(x, e), gpsMeas, gpsParams);
        hMinus = gps_update_model(state_injection(x, -e), gpsMeas, gpsParams);
        Hnum(:,i) = -(hPlus - hMinus) / (2*delta);
    end
    gpsJacErr = max(gpsJacErr, max(abs(H - Hnum), [], 'all'));
end
fprintf('   (reference) gps_update_model max |H - H_fd| = %.2e\n', gpsJacErr);

%% 3. Noise statistics
diffGps_params = diffGps_params_noisy;
clear diffGps_model

nSamples = 20000;
x = random_state();
x(10:13) = [1; 0; 0; 0];        % vertical, body = NED
truthInertial = quaternion_to_rotation(x(10:13))' * baselineBody;
errors = zeros(3, nSamples);
for k = 1:nSamples
    meas = diffGps_model(x, k, frame);
    errors(:,k) = meas.basisVec - truthInertial;
end
sampleStd = std(errors, 0, 2);
reportedStd = sqrt(diag(meas.covariance));
headingErr = atan2(errors(2,:) + truthInertial(2), errors(1,:) + truthInertial(1)) ...
    - atan2(truthInertial(2), truthInertial(1));
headingErr = atan2(sin(headingErr), cos(headingErr));  % wrap to [-pi, pi]
fprintf('3. sample std  = [%.4f %.4f %.4f] m\n', sampleStd);
fprintf('   reported std = [%.4f %.4f %.4f] m\n', reportedStd);
fprintf('   heading error std = %.3f deg on a %.2f m baseline\n', ...
    rad2deg(std(headingErr)), norm(baselineBody));

%% 4. One full update() removes a roll error
updateParams.gateThreshold = 11.345;   % 99%, 3 DOF
updateParams.r1 = frame.r1;
updateParams.r2 = frame.r2;

xNominal = zeros(20,1);
xNominal(10:13) = [1; 0; 0; 0];  % vertical, body = NED
trueError = zeros(19,1);
trueError(10:12) = deg2rad([1; 1; 3]);
xTrue = state_injection(xNominal, trueError);

diffGps_params.noise_std = zeros(3,1);
clear diffGps_model
meas = diffGps_model(xTrue, 1, frame);

[xPost, PPost, accepted] = update(xNominal, P0, @diffGps_update_model, meas, updateParams);

% Rotation about the baseline axis is in the null space of H
[~, H] = diffGps_update_model(xNominal, meas, updateParams);
baselineDir = quaternion_to_rotation(xNominal(10:13))' * baselineBody / norm(baselineBody);
nullErr = norm(H(:,10:12) * baselineDir);

errBefore = rad2deg(attitude_error(xNominal, xTrue));
errAfter = rad2deg(attitude_error(xPost, xTrue));
stdBefore = rad2deg(sqrt(diag(P0(10:12,10:12))));
stdAfter = rad2deg(sqrt(diag(PPost(10:12,10:12))));
fprintf('4. accepted = %d\n', accepted);
fprintf('   attitude error before [x y z] = [%6.3f %6.3f %6.3f] deg\n', errBefore);
fprintf('   attitude error after  [x y z] = [%6.3f %6.3f %6.3f] deg\n', errAfter);
fprintf('   attitude 1-sigma before        = [%6.3f %6.3f %6.3f] deg\n', stdBefore);
fprintf('   attitude 1-sigma after         = [%6.3f %6.3f %6.3f] deg\n', stdAfter);
fprintf('   |H * baseline axis| = %.2e\n', nullErr);

Simulink.data.dictionary.closeAll('-discard');
cd(startDir);

%% Pass/fail
checks = { ...
    'noise-off truth match', truthErr < 1e-12 && lengthErr < 1e-12;
    'Jacobian matches finite difference', jacErr < 1e-6;
    'sign-flip control fails', flippedErr > 1e-2;
    'noise std within 3%', all(abs(sampleStd./reportedStd - 1) < 0.03);
    'update accepted', accepted;
    'roll (z) error removed', abs(errAfter(3)) < 0.05*abs(errBefore(3));
    'pitch (y) error removed', abs(errAfter(2)) < 0.05*abs(errBefore(2));
    'baseline axis in null space of H', nullErr < 1e-12;
    'baseline-axis (x) sigma ~unchanged', abs(stdAfter(1)/stdBefore(1) - 1) < 0.01;
    'heading std ~ noise/baseline', abs(std(headingErr)/(reportedStd(2)/norm(baselineBody)) - 1) < 0.05;
    'roll (z) sigma reduced', stdAfter(3) < 0.1*stdBefore(3)};
fprintf('\n');
for i = 1:size(checks,1)
    fprintf('  [%s] %s\n', ternary(checks{i,2}, 'PASS', 'FAIL'), checks{i,1});
end
assert(all([checks{:,2}]), 'test_diffGps: one or more checks failed');
fprintf('test_diffGps: all checks passed\n');

%% Helpers
function x = random_state()
    q = randn(4,1);
    x = [randn(9,1); q/norm(q); 0.01*randn(7,1)];
end

function phi = attitude_error(xEstimate, xTrue)
% Small-angle rotation vector of C_BI_estimate' * C_BI_true, i.e. the
% dtheta that state_injection would need to reach truth
    E = quaternion_to_rotation(xEstimate(10:13))' * quaternion_to_rotation(xTrue(10:13));
    phi = 0.5 * [E(3,2) - E(2,3); E(1,3) - E(3,1); E(2,1) - E(1,2)];
end

function out = ternary(cond, a, b)
    if cond
        out = a;
    else
        out = b;
    end
end
