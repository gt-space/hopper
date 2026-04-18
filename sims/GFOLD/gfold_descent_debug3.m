%% Debug 3 — Inspect Test 11 solution, then find exact conflict
clear; clc;

m0=127.14; m_dry=119.90; g0=9.80665; Isp=225; alpha=1/(Isp*g0);
T_nom=500*4.44822; Tmin=0.20*T_nom; Tmax=1.10*T_nom;
u_min=Tmin/m0; u_max=Tmax/m_dry;
r0=[0;0;51]; v0=[0;0;0]; rf=[0;0;1]; vf=[0;0;0];
vz_max=7; theta_max=deg2rad(10); g_vec=[0;0;-g0];
N=40; tf=20; dt=tf/N; t=linspace(0,tf,N+1);

%% Run Test 11 and extract solution (sig decoupled from u)
fprintf('--- Test 11 solution extraction ---\n');
cvx_begin quiet
    variable r11(3,N+1)
    variable v11(3,N+1)
    variable u11(3,N+1)
    variable sig11(1,N+1)
    variable z11(1,N+1)
    minimize(0)
    subject to
        r11(:,1)==r0; v11(:,1)==v0; r11(:,N+1)==rf; v11(:,N+1)==vf;
        z11(1)==log(m0);
        for k=1:N
            r11(:,k+1)==r11(:,k)+dt/2*(v11(:,k)+v11(:,k+1));
            v11(:,k+1)==v11(:,k)+dt/2*((u11(:,k)+g_vec)+(u11(:,k+1)+g_vec));
            z11(k+1)==z11(k)-dt/2*alpha*(sig11(k)+sig11(k+1));
        end
        sig11>=u_min; sig11<=u_max;
        z11>=log(m_dry);
        r11(3,:)>=0;
cvx_end
fprintf('Test 11 Status: %s\n', cvx_status);

u11_mag = sqrt(sum(u11.^2,1));
fprintf('\nu11 stats (specific thrust from kinematics):\n');
fprintf('  ||u|| range:  [%.4f, %.4f] m/s^2\n', min(u11_mag), max(u11_mag));
fprintf('  u(3) range:   [%.4f, %.4f] m/s^2\n', min(u11(3,:)), max(u11(3,:)));
fprintf('  sig11 range:  [%.4f, %.4f] m/s^2\n', min(sig11), max(sig11));
fprintf('\nChecking if norm(u) <= sig is satisfiable at each node:\n');
violations = sum(u11_mag > sig11 + 1e-6);
fprintf('  Nodes where ||u|| > sig: %d / %d\n', violations, N+1);
fprintf('  Max excess: %.6f\n', max(u11_mag - sig11));

% Check if u(3) >= sig*cos(theta) is satisfiable
pointing_ok = u11(3,:) >= sig11*cos(theta_max) - 1e-6;
fprintf('\nChecking pointing constraint u(3) >= sig*cos(10deg):\n');
fprintf('  sig*cos(theta) range: [%.4f, %.4f]\n', ...
    min(sig11*cos(theta_max)), max(sig11*cos(theta_max)));
fprintf('  Nodes where u(3) < sig*cos(theta): %d / %d\n', ...
    sum(~pointing_ok), N+1);
fprintf('  Min u(3): %.4f,  Min sig*cos(theta): %.4f\n', ...
    min(u11(3,:)), min(sig11*cos(theta_max)));

% The real question: can we simultaneously have
% norm(u) <= sig AND u(3) >= sig*cos(theta)?
% This requires: u(3)/norm(u) >= cos(theta) i.e. thrust angle <= theta_max
% Check what angle the Test 11 u vectors make
angles = acosd(u11(3,:) ./ max(u11_mag, 1e-6));
fprintf('\nThrust vector angles from vertical:\n');
fprintf('  Range: [%.2f, %.2f] deg\n', min(angles), max(angles));
fprintf('  theta_max = %.1f deg\n', rad2deg(theta_max));
fprintf('  Nodes exceeding theta_max: %d / %d\n', sum(angles > rad2deg(theta_max)+0.01), N+1);

fprintf('\n--- KEY DIAGNOSIS ---\n');
if any(angles > rad2deg(theta_max) + 0.1)
    fprintf('FOUND IT: Test 11 produces thrust vectors > %.1f deg from vertical\n', rad2deg(theta_max));
    fprintf('When we add norm(u)<=sig + pointing, these nodes become infeasible\n');
    fprintf('Solution: RELAX theta_max\n');
    % Find minimum required theta
    min_theta_needed = max(angles);
    fprintf('Minimum theta_max needed: %.2f deg\n', min_theta_needed);
else
    fprintf('Thrust angles are all within theta_max — problem is elsewhere\n');
end