%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test_gmekf_predict.m
% 
% Lowkey made with Claude bc i do not have enough time, but I reviewed it 
% and it looks good - David
% 
% Checks for gmekf_predict:
%   1. Interface: sizes, finite values, unit quaternion, aux contents
%   2. xnp1 matches propagate_state called directly; dt = 0 is identity
%   3. Pnp1 is symmetric and positive semi-definite
%   4. With Qc = 0, Pnp1 = STM*Pn*STM' (using aux.F)
%   5. With Pn = 0 and small dt, each noise-driven block grows like
%      (G*Qc*G')*dt  (checks the G/Qc column layout)
%   6. F is the Jacobian of the real propagation: the error-state STM from
%      finite differences through propagate_state matches I + F*dt + F^2*dt^2/2.
%      Error convention: error = estimate - truth, attitude error is
%      inertial: C_est = (I + [alpha x]) C_true.
%
% Errors if any check fails.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rng(7);
g = 9.81;
J = [0.50 0.02 0.01; 0.02 0.60 0.03; 0.01 0.03 0.80];
qa = 0.2; qt = 0.3; qab = 0.05; qgb = 0.07; qbb = 0.11;
Qc = diag([qa*ones(3,1); qt*ones(3,1); qab*ones(3,1); qgb*ones(3,1); qbb]);

checks = cell(0,2);
add = @(name, ok) evalin('base', sprintf('checks(end+1,:) = {''%s'', %d};', name, ok));

%% 1-5: general behavior over random trials
nTrials = 50;
dt = 0.01;
ifaceOK = true; dt0OK = true; matchOK = true; symOK = true; psdOK = true;
stmOK = true; invalidOK = true;
for k = 1:nTrials
    [xn, Pn, accel, torque] = random_case(g);

    for validity = [true false]
        [xnp1, Pnp1, aux] = gmekf_predict(xn, Pn, accel, torque, validity, J, g, dt, Qc);

        ok = isequal(size(xnp1), [20 1]) && isequal(size(Pnp1), [19 19]) ...
            && all(isfinite(xnp1)) && all(isfinite(Pnp1(:))) ...
            && abs(norm(xnp1(10:13)) - 1) < 1e-9 ...
            && all(isfield(aux, {'A_IB','A_IB05','a_I','w_I','F','G'})) ...
            && isequal(size(aux.F), [19 19]) && isequal(size(aux.G), [19 13]) ...
            && norm(aux.A_IB*aux.A_IB' - eye(3), 'fro') < 1e-9 ...
            && norm(aux.A_IB - quaternion_to_rotation(xn(10:13)), 'fro') < 1e-12;
        ifaceOK = ifaceOK && ok;

        symOK = symOK && norm(Pnp1 - Pnp1', 'fro') < 1e-12 * max(1, norm(Pnp1, 'fro'));
        psdOK = psdOK && min(eig(0.5*(Pnp1 + Pnp1'))) > -1e-10;

        if ~validity
            invalidOK = invalidOK && all(isfinite(aux.a_I));
        end
    end

    % xnp1 equals propagate_state directly
    [xnp1, ~, ~] = gmekf_predict(xn, Pn, accel, torque, true, J, g, dt, Qc);
    auxProp = struct('A_IB', quaternion_to_rotation(xn(10:13)));
    xRef = propagate_state(xn, auxProp, accel, torque, true, J, g, dt);
    matchOK = matchOK && norm(xnp1 - xRef) < 1e-12;

    % dt = 0 is the identity
    [x0, P0, ~] = gmekf_predict(xn, Pn, accel, torque, true, J, g, 0, Qc);
    dt0OK = dt0OK && norm(x0 - xn) < 1e-12 && norm(P0 - Pn, 'fro') < 1e-12;

    % Qc = 0: pure STM propagation
    [~, Pq0, auxq] = gmekf_predict(xn, Pn, accel, torque, true, J, g, dt, zeros(13));
    STM = eye(19) + auxq.F*dt + 0.5*(auxq.F*auxq.F)*dt^2;
    stmOK = stmOK && norm(Pq0 - 0.5*(STM*Pn*STM' + (STM*Pn*STM')'), 'fro') ...
        < 1e-10 * max(1, norm(Pq0, 'fro'));
end
checks(end+1,:) = {'interface (sizes, finite, unit quat, aux fields)', ifaceOK};
checks(end+1,:) = {'xnp1 matches propagate_state', matchOK};
checks(end+1,:) = {'dt = 0 leaves x and P unchanged', dt0OK};
checks(end+1,:) = {'Pnp1 symmetric', symOK};
checks(end+1,:) = {'Pnp1 positive semi-definite', psdOK};
checks(end+1,:) = {'Qc = 0 gives STM*Pn*STM''', stmOK};
checks(end+1,:) = {'invalid-IMU branch gives finite a_I', invalidOK};

%% 5: process noise growth layout (Pn = 0, small dt)
[xn, ~, accel, torque] = random_case(g);
dtS = 1e-4;
[~, Pn1, ~] = gmekf_predict(xn, zeros(19), accel, torque, true, J, g, dtS, Qc);
Ji = inv(J);
exp_vel  = qa*dtS*eye(3);
exp_om   = qt*dtS*(Ji*Ji');
exp_ab   = qab*dtS*eye(3);
exp_gb   = qgb*dtS*eye(3);
exp_baro = qbb*dtS;
relerr = @(A, B) norm(A - B, 'fro') / norm(B, 'fro');
noiseErr = [relerr(Pn1(4:6,4:6),   exp_vel), ...
            relerr(Pn1(7:9,7:9),   exp_om), ...
            relerr(Pn1(13:15,13:15), exp_ab), ...
            relerr(Pn1(16:18,16:18), exp_gb), ...
            abs(Pn1(19,19) - exp_baro)/exp_baro];
fprintf('5. noise growth rel. error [vel, omega, accelBias, gyroBias, baro]: %s\n', ...
    mat2str(noiseErr, 2));
checks(end+1,:) = {'process-noise growth matches G*Qc*G''*dt', all(noiseErr < 0.02)};

%% 6: F against finite differences through propagate_state
dtF = 1e-3;
delta = 1e-5;
stmErr = 0;
for k = 1:10
    [xn, Pn, accel, torque] = random_case(g);
    [~, ~, aux] = gmekf_predict(xn, Pn, accel, torque, true, J, g, dtF, Qc);
    STM = eye(19) + aux.F*dtF + 0.5*(aux.F*aux.F)*dtF^2;

    Phi = zeros(19);
    for i = 1:19
        e = zeros(19,1); e(i) = delta;
        Phi(:,i) = (err_after(xn, e,  accel, torque, J, g, dtF) ...
                  - err_after(xn, -e, accel, torque, J, g, dtF)) / (2*delta);
    end
    stmErr = max(stmErr, max(abs(Phi - STM), [], 'all'));
    if k == 1, PhiFirst = Phi; STMFirst = STM; end
end
fprintf('6. max |STM_fd - STM(F)| = %.1e\n', stmErr);
if stmErr >= 1e-5
    % Per-block report to localize the problem
    blocks = {'pos',1:3; 'vel',4:6; 'omega',7:9; 'att',10:12; ...
              'accBias',13:15; 'gyrBias',16:18; 'baro',19};
    fprintf('   per-block max error (first trial), rows x cols:\n');
    for r = 1:size(blocks,1)
        for c = 1:size(blocks,1)
            d = max(abs(PhiFirst(blocks{r,2},blocks{c,2}) - STMFirst(blocks{r,2},blocks{c,2})), [], 'all');
            if d > 1e-6
                fprintf('   %-8s <- %-8s : %.1e\n', blocks{r,1}, blocks{c,1}, d);
            end
        end
    end
end
checks(end+1,:) = {'F matches finite-difference error-state STM', stmErr < 1e-5};

%% Pass/fail
fprintf('\n');
for i = 1:size(checks,1)
    if checks{i,2}, tag = 'PASS'; else, tag = 'FAIL'; end
    fprintf('  [%s] %s\n', tag, checks{i,1});
end
assert(all([checks{:,2}]), 'test_gmekf_predict: one or more checks failed');
fprintf('test_gmekf_predict: all checks passed\n');

%% Helpers
function [xn, Pn, accel, torque] = random_case(g)
    q = randn(4,1); q = q/norm(q);
    xn = zeros(20,1);
    xn(1:3)   = randn(3,1);
    xn(4:6)   = randn(3,1);
    xn(7:9)   = 0.5*randn(3,1);
    xn(10:13) = q;
    xn(14:20) = 0.01*randn(7,1);
    A = randn(19);
    Pn = 1e-4*(A*A') + 1e-6*eye(19);
    accel  = [0; 0; -g] + 0.3*randn(3,1);   % body-frame specific force
    torque = 0.1*randn(3,1);
end

function e_out = err_after(xn, e, accel, torque, J, g, dt)
% Propagate the estimate and a truth that differs by error e
% (error = estimate - truth), then return the error after dt.
    xTrue = state_injection(xn, e);   % truth = estimate - e
    xEst  = xn;
    est  = prop(xEst,  accel, torque, J, g, dt);
    tru  = prop(xTrue, accel, torque, J, g, dt);

    Ce = quaternion_to_rotation(est(10:13));    % C_IB
    Ct = quaternion_to_rotation(tru(10:13));
    M  = Ce*Ct';                                 % I + [alpha x]
    alpha = 0.5*[M(3,2)-M(2,3); M(1,3)-M(3,1); M(2,1)-M(1,2)];

    e_out = [est(1:3)   - tru(1:3);
             est(4:6)   - tru(4:6);
             est(7:9)   - tru(7:9);
             alpha;
             Ce*(est(14:16) - tru(14:16));       % inertial-aligned accel bias error
             Ce*(est(17:19) - tru(17:19));       % inertial-aligned gyro bias error
             est(20)    - tru(20)];
end

function x = prop(x0, accel, torque, J, g, dt)
    auxP = struct('A_IB', quaternion_to_rotation(x0(10:13)));
    x = propagate_state(x0, auxP, accel, torque, true, J, g, dt);
end