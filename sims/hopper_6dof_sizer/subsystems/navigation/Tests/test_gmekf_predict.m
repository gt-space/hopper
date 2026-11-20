%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test_gmekf_predict.m
%
% Checks the predictor's direct gyro-rate input, torque-free Euler
% propagation when the IMU is invalid, and discrete covariance model.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rng(7);
g = 9.81;
J = [0.50 0.02 0.01; 0.02 0.60 0.03; 0.01 0.03 0.80];
qa = 0.2; qg = 0.3; qab = 0.05; qgb = 0.07; qbb = 0.11;
Qc = diag([qa*ones(3,1); qg*ones(3,1); qab*ones(3,1); ...
           qgb*ones(3,1); qbb]);

checks = cell(0,2);
nTrials = 50;
dt = 0.01;
ifaceOK = true; validRateOK = true; invalidEulerOK = true;
symOK = true; psdOK = true; zeroStepOK = true;
for k = 1:nTrials
    [xn, Pn, accel, gyro] = random_case(g);

    for validity = [true false]
        [xnp1, Pnp1, aux] = gmekf_predict( ...
            xn, Pn, accel, gyro, validity, J, g, dt, Qc);

        ifaceOK = ifaceOK && isequal(size(xnp1), [20 1]) ...
            && isequal(size(Pnp1), [19 19]) ...
            && all(isfinite(xnp1)) && all(isfinite(Pnp1(:))) ...
            && abs(norm(xnp1(10:13)) - 1) < 1e-9 ...
            && all(isfield(aux, {'A_IB','A_IB05','a_I','w_I','F','G'})) ...
            && isequal(size(aux.F), [19 19]) && isequal(size(aux.G), [19 13]) ...
            && norm(aux.A_IB*aux.A_IB' - eye(3), 'fro') < 1e-9;

        symOK = symOK && norm(Pnp1 - Pnp1', 'fro') ...
            < 1e-12 * max(1, norm(Pnp1, 'fro'));
        psdOK = psdOK && min(eig(0.5*(Pnp1 + Pnp1'))) > -1e-10;

        auxProp = struct('A_IB', quaternion_to_rotation(xn(10:13)));
        xRef = propagate_state(xn, auxProp, accel, gyro, validity, J, g, dt);
        ifaceOK = ifaceOK && norm(xnp1 - xRef) < 1e-12;

        if validity
            validRateOK = validRateOK && norm(xnp1(7:9) - ...
                (gyro.angvel - xn(17:19))) < 1e-12;
        else
            omegaDot = -J \ cross(xn(7:9), J*xn(7:9));
            invalidEulerOK = invalidEulerOK && norm(xnp1(7:9) - ...
                (xn(7:9) + omegaDot*dt)) < 1e-12;
        end
    end

    % A zero-length valid step still loads the measured, bias-corrected rate.
    [x0, Pzero, ~] = gmekf_predict( ...
        xn, Pn, accel, gyro, true, J, g, 0, zeros(13));
    zeroStepOK = zeroStepOK && norm(x0(1:6) - xn(1:6)) < 1e-12 ...
        && norm(x0(7:9) - (gyro.angvel - xn(17:19))) < 1e-12 ...
        && norm(x0(10:20) - xn(10:20)) < 1e-12 ...
        && all(isfinite(Pzero(:)));
end
checks(end+1,:) = {'interface, finite output, unit quaternion, aux dimensions', ifaceOK};
checks(end+1,:) = {'valid gyro directly sets rate minus bias', validRateOK};
checks(end+1,:) = {'invalid gyro uses torque-free Euler dynamics', invalidEulerOK};
checks(end+1,:) = {'zero step does not integrate state but loads valid gyro rate', zeroStepOK};
checks(end+1,:) = {'Pnp1 symmetric', symOK};
checks(end+1,:) = {'Pnp1 positive semi-definite', psdOK};

%% Process-noise layout and scaling
[xn, ~, accel, gyro] = random_case(g);
dtS = 1e-4;
[~, Pn1, ~] = gmekf_predict(xn, zeros(19), accel, gyro, true, J, g, dtS, Qc);
noiseErr = [ ...
    relerr(Pn1(4:6,4:6), qa*dtS^2*eye(3)), ...
    relerr(Pn1(7:9,7:9), qg*eye(3)), ...
    relerr(Pn1(13:15,13:15), qab*dtS*eye(3)), ...
    relerr(Pn1(16:18,16:18), qgb*dtS*eye(3)), ...
    abs(Pn1(19,19) - qbb*dtS)/(qbb*dtS)];
fprintf('process-noise relative errors [vel, rate, accelBias, gyroBias, baro]: %s\n', ...
    mat2str(noiseErr, 2));
checks(end+1,:) = {'sample IMU covariance and bias random-walk scaling', all(noiseErr < 1e-10)};

%% Discrete F against finite differences of the state propagation
dtF = 1e-3;
delta = 1e-5;
stmErr = 0;
blocks = {1:3; 4:6; 7:9; 10:12; 13:15; 16:18; 19};
blockErr = zeros(7);
for validity = [true false]
    for k = 1:10
        [xn, Pn, accel, gyro] = random_case(g);
        [~, ~, aux] = gmekf_predict( ...
            xn, Pn, accel, gyro, validity, J, g, dtF, Qc);
        Phi = zeros(19);
        for i = 1:19
            e = zeros(19,1); e(i) = delta;
            Phi(:,i) = (err_after(xn, e, accel, gyro, validity, J, g, dtF) ...
                      - err_after(xn, -e, accel, gyro, validity, J, g, dtF)) ...
                      / (2*delta);
        end
        thisErr = max(abs(Phi - aux.F), [], 'all');
        if thisErr > stmErr
            stmErr = thisErr;
            PhiMax = Phi;
            FMax = aux.F;
            validMax = validity;
        end
        for r = 1:numel(blocks)
            for c = 1:numel(blocks)
                blockErr(r,c) = max(blockErr(r,c), ...
                    max(abs(Phi(blocks{r},blocks{c}) - ...
                    aux.F(blocks{r},blocks{c})), [], 'all'));
            end
        end
    end
end
fprintf('max |finite-difference transition - F| = %.1e\n', stmErr);
if stmErr >= 1e-5
    fprintf('   max block errors (row block x column block):\n');
    disp(blockErr);
    fprintf('   validity = %d, attitude rows vs rate/bias columns:\n', validMax);
    disp([PhiMax(10:12,7:9), FMax(10:12,7:9), ...
          PhiMax(10:12,16:18), FMax(10:12,16:18)]);
    fprintf('   gyro-bias transition:\n');
    disp([PhiMax(16:18,10:12), FMax(16:18,10:12), ...
          PhiMax(16:18,16:18), FMax(16:18,16:18)]);
end
checks(end+1,:) = {'discrete F matches finite-difference state transition', stmErr < 1e-5};

%% Qc = 0 applies the discrete transition exactly
[xn, Pn, accel, gyro] = random_case(g);
[~, Pq0, auxq] = gmekf_predict(xn, Pn, accel, gyro, true, J, g, dt, zeros(13));
Pexpected = auxq.F*Pn*auxq.F';
checks(end+1,:) = {'Qc = 0 gives F*Pn*F''', ...
    norm(Pq0 - 0.5*(Pexpected + Pexpected'), 'fro') ...
    < 1e-12 * max(1, norm(Pq0, 'fro'))};

fprintf('\n');
for i = 1:size(checks,1)
    if checks{i,2}, tag = 'PASS'; else, tag = 'FAIL'; end
    fprintf('  [%s] %s\n', tag, checks{i,1});
end
assert(all([checks{:,2}]), 'test_gmekf_predict: one or more checks failed');
fprintf('test_gmekf_predict: all checks passed\n');

%% Helpers
function [xn, Pn, accel, gyro] = random_case(g)
    q = randn(4,1); q = q/norm(q);
    xn = zeros(20,1);
    xn(1:3)   = randn(3,1);
    xn(4:6)   = randn(3,1);
    xn(7:9)   = 0.5*randn(3,1);
    xn(10:13) = q;
    xn(14:20) = 0.01*randn(7,1);
    A = randn(19);
    Pn = 1e-4*(A*A') + 1e-6*eye(19);
    accel = [0; 0; -g] + 0.3*randn(3,1);
    gyro = struct('angvel', 0.5*randn(3,1));
end

function e_out = err_after(xn, e, accel, gyro, validity, J, g, dt)
    xTrue = state_injection(xn, e);
    est = prop(xn, accel, gyro, validity, J, g, dt);
    tru = prop(xTrue, accel, gyro, validity, J, g, dt);

    Ce = quaternion_to_rotation(est(10:13));
    Ct = quaternion_to_rotation(tru(10:13));
    M = Ce*Ct';
    alpha = 0.5*[M(3,2)-M(2,3); M(1,3)-M(3,1); M(2,1)-M(1,2)];

    e_out = [est(1:3) - tru(1:3);
             est(4:6) - tru(4:6);
             est(7:9) - tru(7:9);
             alpha;
             Ce*(est(14:16) - tru(14:16));
             Ce*(est(17:19) - tru(17:19));
             est(20) - tru(20)];
end

function x = prop(x0, accel, gyro, validity, J, g, dt)
    auxP = struct('A_IB', quaternion_to_rotation(x0(10:13)));
    x = propagate_state(x0, auxP, accel, gyro, validity, J, g, dt);
end

function err = relerr(A, B)
    err = norm(A - B, 'fro') / max(norm(B, 'fro'), eps);
end
