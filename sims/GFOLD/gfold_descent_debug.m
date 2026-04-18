% Minimal descent debug — strip everything to bare bones
clear; clc;

% Params
m0    = 127.14; m_dry = 119.90; g0 = 9.80665;
T_nom = 500*4.44822; Tmin = 0.20*T_nom; Tmax = 1.10*T_nom;
Isp   = 225; alpha = 1/(Isp*g0);
g_vec = [0;0;-g0];
r0 = [0;0;51]; v0 = [0;0;0];
rf = [0;0;1];  vf = [0;0;0];
N = 40; tf = 20; dt = tf/N;
t = linspace(0,tf,N+1);

u_min = Tmin/m0;
u_max = Tmax/m_dry;

fprintf('u_min=%.3f u_max=%.3f g=%.3f\n', u_min, u_max, g0);
fprintf('Net accel range: [%.3f, %.3f] m/s^2\n', u_min-g0, u_max-g0);

%% Test 1: No mass dynamics, no thrust bounds, just kinematics
fprintf('\n--- Test 1: Pure kinematics only ---\n');
cvx_begin quiet
    variable r1(3,N+1)
    variable v1(3,N+1)
    variable u1(3,N+1)
    minimize(0)
    subject to
        r1(:,1) == r0; v1(:,1) == v0;
        r1(:,N+1) == rf; v1(:,N+1) == vf;
        for k = 1:N
            r1(:,k+1) == r1(:,k) + dt/2*(v1(:,k)+v1(:,k+1));
            v1(:,k+1) == v1(:,k) + dt/2*((u1(:,k)+g_vec)+(u1(:,k+1)+g_vec));
        end
        r1(3,:) >= 0;
cvx_end
fprintf('Status: %s\n', cvx_status);

%% Test 2: Add thrust magnitude bounds
fprintf('\n--- Test 2: + thrust bounds ---\n');
cvx_begin quiet
    variable r2(3,N+1)
    variable v2(3,N+1)
    variable u2(3,N+1)
    variable sig2(1,N+1)
    minimize(0)
    subject to
        r2(:,1) == r0; v2(:,1) == v0;
        r2(:,N+1) == rf; v2(:,N+1) == vf;
        for k = 1:N
            r2(:,k+1) == r2(:,k) + dt/2*(v2(:,k)+v2(:,k+1));
            v2(:,k+1) == v2(:,k) + dt/2*((u2(:,k)+g_vec)+(u2(:,k+1)+g_vec));
        end
        for k = 1:N+1
            norm(u2(:,k),2) <= sig2(k);
            sig2(k) >= u_min;
            sig2(k) <= u_max;
        end
        r2(3,:) >= 0;
        sig2 >= 0;
cvx_end
fprintf('Status: %s\n', cvx_status);

%% Test 3: Add pointing constraint
fprintf('\n--- Test 3: + pointing constraint (theta=10 deg) ---\n');
theta_max = deg2rad(10);
cvx_begin quiet
    variable r3(3,N+1)
    variable v3(3,N+1)
    variable u3(3,N+1)
    variable sig3(1,N+1)
    minimize(0)
    subject to
        r3(:,1) == r0; v3(:,1) == v0;
        r3(:,N+1) == rf; v3(:,N+1) == vf;
        for k = 1:N
            r3(:,k+1) == r3(:,k) + dt/2*(v3(:,k)+v3(:,k+1));
            v3(:,k+1) == v3(:,k) + dt/2*((u3(:,k)+g_vec)+(u3(:,k+1)+g_vec));
        end
        for k = 1:N+1
            norm(u3(:,k),2) <= sig3(k);
            sig3(k) >= u_min;
            sig3(k) <= u_max;
            u3(3,k) >= sig3(k)*cos(theta_max);
        end
        r3(3,:) >= 0;
        sig3 >= 0;
cvx_end
fprintf('Status: %s\n', cvx_status);

%% Test 4: Add vz rate limit
fprintf('\n--- Test 4: + vz rate limit (7 m/s) ---\n');
vz_max = 7;
cvx_begin quiet
    variable r4(3,N+1)
    variable v4(3,N+1)
    variable u4(3,N+1)
    variable sig4(1,N+1)
    minimize(0)
    subject to
        r4(:,1) == r0; v4(:,1) == v0;
        r4(:,N+1) == rf; v4(:,N+1) == vf;
        for k = 1:N
            r4(:,k+1) == r4(:,k) + dt/2*(v4(:,k)+v4(:,k+1));
            v4(:,k+1) == v4(:,k) + dt/2*((u4(:,k)+g_vec)+(u4(:,k+1)+g_vec));
        end
        for k = 1:N+1
            norm(u4(:,k),2) <= sig4(k);
            sig4(k) >= u_min;
            sig4(k) <= u_max;
            u4(3,k) >= sig4(k)*cos(theta_max);
            v4(3,k) >= -vz_max;
            v4(3,k) <= vz_max;
        end
        r4(3,:) >= 0;
        sig4 >= 0;
cvx_end
fprintf('Status: %s\n', cvx_status);

%% Test 5: + mass dynamics, no objective
fprintf('\n--- Test 5: + mass dynamics (z variable), feasibility only ---\n');
vz_max = 7; theta_max = deg2rad(10);
alpha = 1/(225*9.80665);
cvx_begin quiet
    variable r5(3,N+1)
    variable v5(3,N+1)
    variable u5(3,N+1)
    variable sig5(1,N+1)
    variable z5(1,N+1)
    minimize(0)
    subject to
        r5(:,1) == r0; v5(:,1) == v0;
        r5(:,N+1) == rf; v5(:,N+1) == vf;
        z5(1) == log(m0);
        for k = 1:N
            r5(:,k+1) == r5(:,k) + dt/2*(v5(:,k)+v5(:,k+1));
            v5(:,k+1) == v5(:,k) + dt/2*((u5(:,k)+g_vec)+(u5(:,k+1)+g_vec));
            z5(k+1) == z5(k) - dt/2*alpha*(sig5(k)+sig5(k+1));
        end
        for k = 1:N+1
            norm(u5(:,k),2) <= sig5(k);
            sig5(k) >= u_min;
            sig5(k) <= u_max;
            u5(3,k) >= sig5(k)*cos(theta_max);
            v5(3,k) >= -vz_max;
            v5(3,k) <= vz_max;
        end
        for k = 1:N+1
            z5(k) >= log(m_dry);
            z5(k) <= log(m0);
        end
        r5(3,:) >= 0;
        sig5 >= 0;
cvx_end
fprintf('Status: %s\n', cvx_status);

%% Test 6: + maximize final mass (full objective)
fprintf('\n--- Test 6: + maximize z(N+1) objective ---\n');
cvx_begin quiet
    variable r6(3,N+1)
    variable v6(3,N+1)
    variable u6(3,N+1)
    variable sig6(1,N+1)
    variable z6(1,N+1)
    maximize(z6(N+1))
    subject to
        r6(:,1) == r0; v6(:,1) == v0;
        r6(:,N+1) == rf; v6(:,N+1) == vf;
        z6(1) == log(m0);
        for k = 1:N
            r6(:,k+1) == r6(:,k) + dt/2*(v6(:,k)+v6(:,k+1));
            v6(:,k+1) == v6(:,k) + dt/2*((u6(:,k)+g_vec)+(u6(:,k+1)+g_vec));
            z6(k+1) == z6(k) - dt/2*alpha*(sig6(k)+sig6(k+1));
        end
        for k = 1:N+1
            norm(u6(:,k),2) <= sig6(k);
            sig6(k) >= u_min;
            sig6(k) <= u_max;
            u6(3,k) >= sig6(k)*cos(theta_max);
            v6(3,k) >= -vz_max;
            v6(3,k) <= vz_max;
        end
        for k = 1:N+1
            z6(k) >= log(m_dry);
            z6(k) <= log(m0);
        end
        r6(3,:) >= 0;
        sig6 >= 0;
cvx_end
fprintf('Status: %s\n', cvx_status);
if strcmp(cvx_status,'Solved') || strcmp(cvx_status,'Inaccurate/Solved')
    mass_final = exp(z6(end));
    fuel_used  = m0 - mass_final;
    fprintf('Final mass: %.2f kg, Fuel used: %.2f kg\n', mass_final, fuel_used);
end

%% Test 7: Mass dynamics — drop z <= log(m0) upper bound (redundant)
fprintf('\n--- Test 7: mass dynamics, drop z upper bound ---\n');
vz_max = 7; theta_max = deg2rad(10);
alpha = 1/(225*9.80665);
cvx_begin quiet
    variable r7(3,N+1)
    variable v7(3,N+1)
    variable u7(3,N+1)
    variable sig7(1,N+1)
    variable z7(1,N+1)
    minimize(0)
    subject to
        r7(:,1) == r0; v7(:,1) == v0;
        r7(:,N+1) == rf; v7(:,N+1) == vf;
        z7(1) == log(m0);
        for k = 1:N
            r7(:,k+1) == r7(:,k) + dt/2*(v7(:,k)+v7(:,k+1));
            v7(:,k+1) == v7(:,k) + dt/2*((u7(:,k)+g_vec)+(u7(:,k+1)+g_vec));
            z7(k+1) == z7(k) - dt/2*alpha*(sig7(k)+sig7(k+1));
        end
        for k = 1:N+1
            norm(u7(:,k),2) <= sig7(k);
            sig7(k) >= u_min;
            sig7(k) <= u_max;
            u7(3,k) >= sig7(k)*cos(theta_max);
            v7(3,k) >= -vz_max;
            v7(3,k) <= vz_max;
        end
        z7 >= log(m_dry);   % only lower bound — upper bound is redundant (mass only decreases)
        r7(3,:) >= 0;
        sig7 >= 0;
cvx_end
fprintf('Status: %s\n', cvx_status);

%% Test 8: Same as 7 but with maximize objective
fprintf('\n--- Test 8: drop z upper bound + maximize ---\n');
cvx_begin quiet
    variable r8(3,N+1)
    variable v8(3,N+1)
    variable u8(3,N+1)
    variable sig8(1,N+1)
    variable z8(1,N+1)
    maximize(z8(N+1))
    subject to
        r8(:,1) == r0; v8(:,1) == v0;
        r8(:,N+1) == rf; v8(:,N+1) == vf;
        z8(1) == log(m0);
        for k = 1:N
            r8(:,k+1) == r8(:,k) + dt/2*(v8(:,k)+v8(:,k+1));
            v8(:,k+1) == v8(:,k) + dt/2*((u8(:,k)+g_vec)+(u8(:,k+1)+g_vec));
            z8(k+1) == z8(k) - dt/2*alpha*(sig8(k)+sig8(k+1));
        end
        for k = 1:N+1
            norm(u8(:,k),2) <= sig8(k);
            sig8(k) >= u_min;
            sig8(k) <= u_max;
            u8(3,k) >= sig8(k)*cos(theta_max);
            v8(3,k) >= -vz_max;
            v8(3,k) <= vz_max;
        end
        z8 >= log(m_dry);
        r8(3,:) >= 0;
        sig8 >= 0;
cvx_end
fprintf('Status: %s\n', cvx_status);
if strcmp(cvx_status,'Solved') || strcmp(cvx_status,'Inaccurate/Solved')
    mass_final = exp(z8(end));
    fprintf('Final mass: %.2f kg, Fuel used: %.2f kg\n', mass_final, m0-mass_final);
end

%% Test 9: Mass dynamics, no sig lower bound
fprintf('\n--- Test 9: mass dynamics, no sig lower bound ---\n');
vz_max = 7; theta_max = deg2rad(10);
alpha = 1/(225*9.80665);
cvx_begin quiet
    variable r9(3,N+1)
    variable v9(3,N+1)
    variable u9(3,N+1)
    variable sig9(1,N+1)
    variable z9(1,N+1)
    minimize(0)
    subject to
        r9(:,1) == r0; v9(:,1) == v0;
        r9(:,N+1) == rf; v9(:,N+1) == vf;
        z9(1) == log(m0);
        for k = 1:N
            r9(:,k+1) == r9(:,k) + dt/2*(v9(:,k)+v9(:,k+1));
            v9(:,k+1) == v9(:,k) + dt/2*((u9(:,k)+g_vec)+(u9(:,k+1)+g_vec));
            z9(k+1) == z9(k) - dt/2*alpha*(sig9(k)+sig9(k+1));
        end
        for k = 1:N+1
            norm(u9(:,k),2) <= sig9(k);
            % NO sig lower bound
            sig9(k) <= u_max;
            u9(3,k) >= sig9(k)*cos(theta_max);
            v9(3,k) >= -vz_max;
            v9(3,k) <= vz_max;
        end
        z9 >= log(m_dry);
        r9(3,:) >= 0;
        sig9 >= 0;
cvx_end
fprintf('Status: %s\n', cvx_status);

%% Test 10: Mass dynamics only, bare minimum — just z equality chain + bounds
fprintf('\n--- Test 10: z equality chain only, no kinematics ---\n');
cvx_begin quiet
    variable sig10(1,N+1)
    variable z10(1,N+1)
    minimize(0)
    subject to
        z10(1) == log(m0);
        for k = 1:N
            z10(k+1) == z10(k) - dt/2*alpha*(sig10(k)+sig10(k+1));
        end
        sig10 >= u_min;
        sig10 <= u_max;
        z10 >= log(m_dry);
cvx_end
fprintf('Status: %s\n', cvx_status);
if strcmp(cvx_status,'Solved') || strcmp(cvx_status,'Inaccurate/Solved')
    fprintf('z range: [%.4f, %.4f], log(m_dry)=%.4f\n', min(z10), max(z10), log(m_dry));
    fprintf('Mass range: [%.2f, %.2f] kg\n', exp(min(z10)), exp(max(z10)));
end

%% Test 11: z chain + kinematics, no sig coupling to u
fprintf('\n--- Test 11: z chain + kinematics, sig decoupled from u ---\n');
cvx_begin quiet
    variable r11(3,N+1)
    variable v11(3,N+1)
    variable u11(3,N+1)
    variable sig11(1,N+1)
    variable z11(1,N+1)
    minimize(0)
    subject to
        r11(:,1) == r0; v11(:,1) == v0;
        r11(:,N+1) == rf; v11(:,N+1) == vf;
        z11(1) == log(m0);
        for k = 1:N
            r11(:,k+1) == r11(:,k) + dt/2*(v11(:,k)+v11(:,k+1));
            v11(:,k+1) == v11(:,k) + dt/2*((u11(:,k)+g_vec)+(u11(:,k+1)+g_vec));
            z11(k+1) == z11(k) - dt/2*alpha*(sig11(k)+sig11(k+1));
        end
        % sig and u are independent — no norm(u)<=sig coupling
        sig11 >= u_min;
        sig11 <= u_max;
        z11 >= log(m_dry);
        r11(3,:) >= 0;
cvx_end
fprintf('Status: %s\n', cvx_status);

%% Test 12: Split sig (mass) from T_mag (thrust SOC) — key fix
fprintf('\n--- Test 12: split sig/T_mag, full problem ---\n');
vz_max = 7; theta_max = deg2rad(10);
alpha = 1/(225*9.80665);
cvx_begin quiet
    variable r12(3,N+1)
    variable v12(3,N+1)
    variable u12(3,N+1)
    variable sig12(1,N+1)   % drives mass depletion only
    variable z12(1,N+1)
    minimize(0)
    subject to
        r12(:,1) == r0; v12(:,1) == v0;
        r12(:,N+1) == rf; v12(:,N+1) == vf;
        z12(1) == log(m0);
        for k = 1:N
            r12(:,k+1) == r12(:,k) + dt/2*(v12(:,k)+v12(:,k+1));
            v12(:,k+1) == v12(:,k) + dt/2*((u12(:,k)+g_vec)+(u12(:,k+1)+g_vec));
            z12(k+1)   == z12(k) - dt/2*alpha*(sig12(k)+sig12(k+1));
        end
        for k = 1:N+1
            % SOC: ||u|| <= sig (lossless convexification)
            norm(u12(:,k),2) <= sig12(k);
            % Thrust bounds via sig
            sig12(k) >= u_min;
            sig12(k) <= u_max;
            % Pointing: thrust must be upward within theta_max
            u12(3,k) >= sig12(k) * cos(theta_max);
            % Velocity limits
            v12(3,k) >= -vz_max;
            v12(3,k) <= vz_max;
        end
        z12 >= log(m_dry);
        r12(3,:) >= 0;
        sig12 >= 0;
cvx_end
fprintf('Status: %s\n', cvx_status);

%% Test 13: Same as 12 with maximize
fprintf('\n--- Test 13: split sig/T_mag + maximize ---\n');
cvx_begin quiet
    variable r13(3,N+1)
    variable v13(3,N+1)
    variable u13(3,N+1)
    variable sig13(1,N+1)
    variable z13(1,N+1)
    maximize(z13(N+1))
    subject to
        r13(:,1) == r0; v13(:,1) == v0;
        r13(:,N+1) == rf; v13(:,N+1) == vf;
        z13(1) == log(m0);
        for k = 1:N
            r13(:,k+1) == r13(:,k) + dt/2*(v13(:,k)+v13(:,k+1));
            v13(:,k+1) == v13(:,k) + dt/2*((u13(:,k)+g_vec)+(u13(:,k+1)+g_vec));
            z13(k+1)   == z13(k) - dt/2*alpha*(sig13(k)+sig13(k+1));
        end
        for k = 1:N+1
            norm(u13(:,k),2) <= sig13(k);
            sig13(k) >= u_min;
            sig13(k) <= u_max;
            u13(3,k) >= sig13(k)*cos(theta_max);
            v13(3,k) >= -vz_max;
            v13(3,k) <= vz_max;
        end
        z13 >= log(m_dry);
        r13(3,:) >= 0;
        sig13 >= 0;
cvx_end
fprintf('Status: %s\n', cvx_status);
if strcmp(cvx_status,'Solved') || strcmp(cvx_status,'Inaccurate/Solved')
    mass_f = exp(z13(end));
    fprintf('Final mass: %.2f kg, Fuel used: %.2f kg\n', mass_f, m0-mass_f);
    fprintf('Max |vz|: %.3f m/s\n', max(abs(v13(3,:))));
end

%% Test 14: Relax pointing constraint — maybe cos(theta)*sig conflicts with norm(u)<=sig+z chain
fprintf('\n--- Test 14: no pointing constraint, full z chain ---\n');
vz_max = 7;
alpha = 1/(225*9.80665);
cvx_begin quiet
    variable r14(3,N+1)
    variable v14(3,N+1)
    variable u14(3,N+1)
    variable sig14(1,N+1)
    variable z14(1,N+1)
    minimize(0)
    subject to
        r14(:,1) == r0; v14(:,1) == v0;
        r14(:,N+1) == rf; v14(:,N+1) == vf;
        z14(1) == log(m0);
        for k = 1:N
            r14(:,k+1) == r14(:,k) + dt/2*(v14(:,k)+v14(:,k+1));
            v14(:,k+1) == v14(:,k) + dt/2*((u14(:,k)+g_vec)+(u14(:,k+1)+g_vec));
            z14(k+1)   == z14(k) - dt/2*alpha*(sig14(k)+sig14(k+1));
        end
        for k = 1:N+1
            norm(u14(:,k),2) <= sig14(k);
            sig14(k) >= u_min;
            sig14(k) <= u_max;
            % NO pointing constraint
            v14(3,k) >= -vz_max;
            v14(3,k) <= vz_max;
        end
        z14 >= log(m_dry);
        r14(3,:) >= 0;
        sig14 >= 0;
cvx_end
fprintf('Status: %s\n', cvx_status);

%% Test 15: No vz limit, full z chain + pointing
fprintf('\n--- Test 15: no vz limit, full z chain + pointing ---\n');
theta_max = deg2rad(10);
cvx_begin quiet
    variable r15(3,N+1)
    variable v15(3,N+1)
    variable u15(3,N+1)
    variable sig15(1,N+1)
    variable z15(1,N+1)
    minimize(0)
    subject to
        r15(:,1) == r0; v15(:,1) == v0;
        r15(:,N+1) == rf; v15(:,N+1) == vf;
        z15(1) == log(m0);
        for k = 1:N
            r15(:,k+1) == r15(:,k) + dt/2*(v15(:,k)+v15(:,k+1));
            v15(:,k+1) == v15(:,k) + dt/2*((u15(:,k)+g_vec)+(u15(:,k+1)+g_vec));
            z15(k+1)   == z15(k) - dt/2*alpha*(sig15(k)+sig15(k+1));
        end
        for k = 1:N+1
            norm(u15(:,k),2) <= sig15(k);
            sig15(k) >= u_min;
            sig15(k) <= u_max;
            u15(3,k) >= sig15(k)*cos(theta_max);
            % NO vz limit
        end
        z15 >= log(m_dry);
        r15(3,:) >= 0;
        sig15 >= 0;
cvx_end
fprintf('Status: %s\n', cvx_status);

%% Test 16: Replace log-mass z with direct mass variable m
% Mass dynamics: m(k+1) = m(k) - dt/2*alpha*(||T(k)|| + ||T(k+1)||)
% Thrust vector T = u*m, so ||T|| = sig*m (where sig = ||u||)
% But that's nonlinear. Instead track T directly as decision variable.
% m(k+1) = m(k) - dt/2*(alpha_T(k) + alpha_T(k+1))
% where alpha_T = ||T||/Isp/g0, and ||T|| is a separate scalar variable.
fprintf('\n--- Test 16: direct mass variable, T as decision ---\n');
vz_max = 7; theta_max = deg2rad(10);
Isp_val = 225; g0_val = 9.80665;
Tmin_N = 0.20 * 500*4.44822;  % actual thrust in Newtons
Tmax_N = 1.10 * 500*4.44822;

cvx_begin quiet
    variable r16(3,N+1)   % position
    variable v16(3,N+1)   % velocity  
    variable T16(3,N+1)   % thrust vector [N]
    variable Tm16(1,N+1)  % thrust magnitude [N] (slack)
    variable m16(1,N+1)   % mass [kg]
    minimize(0)
    subject to
        r16(:,1) == r0; v16(:,1) == v0; m16(1) == m0;
        r16(:,N+1) == rf; v16(:,N+1) == vf;
        for k = 1:N
            % Kinematics: acceleration = T/m + g
            % Use midpoint mass for division — approximate with m at node k
            % T/m linearized: use fixed m0 as scale (conservative)
            r16(:,k+1) == r16(:,k) + dt/2*(v16(:,k)+v16(:,k+1));
            v16(:,k+1) == v16(:,k) + dt/2*( ...
                (T16(:,k)/m0 + g_vec) + (T16(:,k+1)/m0 + g_vec));
            % Mass depletion: mdot = -||T||/(Isp*g0)
            m16(k+1) == m16(k) - dt/(2*Isp_val*g0_val)*(Tm16(k)+Tm16(k+1));
        end
        for k = 1:N+1
            norm(T16(:,k),2) <= Tm16(k);   % SOC
            Tm16(k) >= Tmin_N;
            Tm16(k) <= Tmax_N;
            T16(3,k) >= Tm16(k)*cos(theta_max);  % pointing
            v16(3,k) >= -vz_max;
            v16(3,k) <= vz_max;
            m16(k)   >= m_dry;
            m16(k)   <= m0;
        end
        r16(3,:) >= 0;
        Tm16 >= 0;
cvx_end
fprintf('Status: %s\n', cvx_status);
if strcmp(cvx_status,'Solved') || strcmp(cvx_status,'Inaccurate/Solved')
    fprintf('Final mass: %.2f kg, Fuel used: %.2f kg\n', m16(end), m0-m16(end));
    fprintf('Max |vz|: %.3f m/s\n', max(abs(v16(3,:))));
end

%% Test 17: Same as Test 5 but force SeDuMi solver
fprintf('\n--- Test 17: Test 5 formulation with SeDuMi solver ---\n');
vz_max = 7; theta_max = deg2rad(10);
alpha = 1/(225*9.80665);
cvx_begin quiet
    cvx_solver sedumi
    variable r17(3,N+1)
    variable v17(3,N+1)
    variable u17(3,N+1)
    variable sig17(1,N+1)
    variable z17(1,N+1)
    minimize(0)
    subject to
        r17(:,1) == r0; v17(:,1) == v0;
        r17(:,N+1) == rf; v17(:,N+1) == vf;
        z17(1) == log(m0);
        for k = 1:N
            r17(:,k+1) == r17(:,k) + dt/2*(v17(:,k)+v17(:,k+1));
            v17(:,k+1) == v17(:,k) + dt/2*((u17(:,k)+g_vec)+(u17(:,k+1)+g_vec));
            z17(k+1)   == z17(k) - dt/2*alpha*(sig17(k)+sig17(k+1));
        end
        for k = 1:N+1
            norm(u17(:,k),2) <= sig17(k);
            sig17(k) >= u_min;
            sig17(k) <= u_max;
            u17(3,k) >= sig17(k)*cos(theta_max);
            v17(3,k) >= -vz_max;
            v17(3,k) <= vz_max;
        end
        z17 >= log(m_dry);
        r17(3,:) >= 0;
        sig17 >= 0;
cvx_end
fprintf('Status: %s\n', cvx_status);
if strcmp(cvx_status,'Solved') || strcmp(cvx_status,'Inaccurate/Solved')
    mass_f = exp(z17(end));
    fprintf('Final mass: %.2f kg, Fuel used: %.2f kg\n', mass_f, m0-mass_f);
end

%% Test 18: SeDuMi + maximize
fprintf('\n--- Test 18: maximize z(N+1) with SeDuMi ---\n');
cvx_begin quiet
    cvx_solver sedumi
    variable r18(3,N+1)
    variable v18(3,N+1)
    variable u18(3,N+1)
    variable sig18(1,N+1)
    variable z18(1,N+1)
    maximize(z18(N+1))
    subject to
        r18(:,1) == r0; v18(:,1) == v0;
        r18(:,N+1) == rf; v18(:,N+1) == vf;
        z18(1) == log(m0);
        for k = 1:N
            r18(:,k+1) == r18(:,k) + dt/2*(v18(:,k)+v18(:,k+1));
            v18(:,k+1) == v18(:,k) + dt/2*((u18(:,k)+g_vec)+(u18(:,k+1)+g_vec));
            z18(k+1)   == z18(k) - dt/2*alpha*(sig18(k)+sig18(k+1));
        end
        for k = 1:N+1
            norm(u18(:,k),2) <= sig18(k);
            sig18(k) >= u_min;
            sig18(k) <= u_max;
            u18(3,k) >= sig18(k)*cos(theta_max);
            v18(3,k) >= -vz_max;
            v18(3,k) <= vz_max;
        end
        z18 >= log(m_dry);
        r18(3,:) >= 0;
        sig18 >= 0;
cvx_end
fprintf('Status: %s\n', cvx_status);
if strcmp(cvx_status,'Solved') || strcmp(cvx_status,'Inaccurate/Solved')
    mass_f = exp(z18(end));
    fprintf('Final mass: %.2f kg, Fuel used: %.2f kg\n', mass_f, m0-mass_f);
    fprintf('Max |vz|: %.3f m/s\n', max(abs(v18(3,:))));
end