%% GFOLD Descent Debug 2 — Proper variable scaling
clear; clc;

%% Physical params
m0    = 127.14;  m_dry = 119.90;
g0    = 9.80665; Isp   = 225;
alpha = 1/(Isp*g0);
T_nom = 500*4.44822;
Tmin  = 0.20*T_nom;  Tmax = 1.10*T_nom;
r0 = [0;0;51]; v0 = [0;0;0];
rf = [0;0;1];  vf = [0;0;0];
vz_max = 7;  theta_max = deg2rad(10);
g_vec = [0;0;-g0];
N = 40;  tf = 20;  dt = tf/N;
t = linspace(0,tf,N+1);

%% Scale factors — normalize all variables to O(1)
Lr = 51;          % position scale [m]
Lv = sqrt(2*g0*51); % velocity scale [m/s] ~ free-fall velocity ~ 31.6
Lt = Lr/Lv;       % time scale [s]
Lm = m0;          % mass scale [kg]
La = g0;          % acceleration scale [m/s^2]

fprintf('Scale factors: Lr=%.1f m, Lv=%.2f m/s, Lt=%.2f s, La=%.3f m/s2\n', Lr,Lv,Lt,La);
fprintf('Nondim tf = %.3f\n', tf/Lt);

% Nondimensional quantities
r0n   = r0/Lr;    rfn   = rf/Lr;
v0n   = v0/Lv;    vfn   = vf/Lv;
g_nd  = g_vec/La;
tfn   = tf/Lt;    dtn   = tfn/N;
u_min_nd = (Tmin/m0)/La;
u_max_nd = (Tmax/m_dry)/La;
alpha_nd = alpha * La * Lt;   % nondim mass flow
mdry_nd  = m_dry/Lm;

fprintf('u_min_nd=%.4f  u_max_nd=%.4f\n', u_min_nd, u_max_nd);
fprintf('alpha_nd=%.6f\n', alpha_nd);

%% Test A: Full problem, scaled, SeDuMi
fprintf('\n--- Test A: full problem, scaled, SeDuMi ---\n');
cvx_begin quiet
    cvx_solver sedumi
    variable r(3,N+1)
    variable v(3,N+1)
    variable u(3,N+1)
    variable sig(1,N+1)
    variable z(1,N+1)
    maximize(z(N+1))
    subject to
        r(:,1) == r0n;  v(:,1) == v0n;  z(1) == log(1.0);
        r(:,N+1) == rfn; v(:,N+1) == vfn;
        for k = 1:N
            r(:,k+1) == r(:,k) + dtn/2*(v(:,k)+v(:,k+1));
            v(:,k+1) == v(:,k) + dtn/2*((u(:,k)+g_nd)+(u(:,k+1)+g_nd));
            z(k+1)   == z(k) - dtn/2*alpha_nd*(sig(k)+sig(k+1));
        end
        for k = 1:N+1
            norm(u(:,k),2) <= sig(k);
            sig(k) >= u_min_nd;
            sig(k) <= u_max_nd;
            u(3,k) >= sig(k)*cos(theta_max);
            v(3,k) >= -vz_max/Lv;
            v(3,k) <= vz_max/Lv;
        end
        z >= log(mdry_nd);
        r(3,:) >= 0;
        sig >= 0;
cvx_end
fprintf('Status: %s\n', cvx_status);
if strcmp(cvx_status,'Solved') || strcmp(cvx_status,'Inaccurate/Solved')
    mass_f = m0*exp(z(end));
    fprintf('Final mass: %.2f kg, Fuel used: %.2f kg\n', mass_f, m0-mass_f);
    fprintf('Max |vz|: %.3f m/s\n', max(abs(v(3,:)))*Lv);
    fprintf('Final pos: [%.3f, %.3f, %.3f] m\n', r(1,end)*Lr, r(2,end)*Lr, r(3,end)*Lr);
end

%% Test B: Same but SDPT3
fprintf('\n--- Test B: full problem, scaled, SDPT3 ---\n');
cvx_begin quiet
    cvx_solver sdpt3
    variable r(3,N+1)
    variable v(3,N+1)
    variable u(3,N+1)
    variable sig(1,N+1)
    variable z(1,N+1)
    maximize(z(N+1))
    subject to
        r(:,1) == r0n;  v(:,1) == v0n;  z(1) == log(1.0);
        r(:,N+1) == rfn; v(:,N+1) == vfn;
        for k = 1:N
            r(:,k+1) == r(:,k) + dtn/2*(v(:,k)+v(:,k+1));
            v(:,k+1) == v(:,k) + dtn/2*((u(:,k)+g_nd)+(u(:,k+1)+g_nd));
            z(k+1)   == z(k) - dtn/2*alpha_nd*(sig(k)+sig(k+1));
        end
        for k = 1:N+1
            norm(u(:,k),2) <= sig(k);
            sig(k) >= u_min_nd;
            sig(k) <= u_max_nd;
            u(3,k) >= sig(k)*cos(theta_max);
            v(3,k) >= -vz_max/Lv;
            v(3,k) <= vz_max/Lv;
        end
        z >= log(mdry_nd);
        r(3,:) >= 0;
        sig >= 0;
cvx_end
fprintf('Status: %s\n', cvx_status);
if strcmp(cvx_status,'Solved') || strcmp(cvx_status,'Inaccurate/Solved')
    mass_f = m0*exp(z(end));
    fprintf('Final mass: %.2f kg, Fuel used: %.2f kg\n', mass_f, m0-mass_f);
    fprintf('Max |vz|: %.3f m/s\n', max(abs(v(3,:)))*Lv);
end