%% Debug 5 — Absolute minimum: 1D vertical descent
clear; clc;

g0=9.80665; m0=127.14; m_dry=119.90; Isp=225; alpha=1/(Isp*g0);
Tmax=1.10*500*4.44822; u_max=Tmax/m_dry;
N=40; tf=20; dt=tf/N; t=linspace(0,tf,N+1);
g_vec=[0;0;-g0];

%% Test A: 1D, no mass, no bounds — just fall from 51 to 1, rest to rest
fprintf('--- Test A: 1D kinematics only, z-axis ---\n');
cvx_begin quiet
    variable z1(1,N+1)   % altitude
    variable vz(1,N+1)   % vertical velocity
    variable uz(1,N+1)   % vertical thrust accel
    minimize(0)
    subject to
        z1(1)==51; vz(1)==0;
        z1(N+1)==1; vz(N+1)==0;
        for k=1:N
            z1(k+1)==z1(k)+dt/2*(vz(k)+vz(k+1));
            vz(k+1)==vz(k)+dt/2*((uz(k)-g0)+(uz(k+1)-g0));
        end
        z1>=0;
cvx_end
fprintf('Status: %s\n', cvx_status);
if strcmp(cvx_status,'Solved')||strcmp(cvx_status,'Inaccurate/Solved')
    fprintf('uz range: [%.3f, %.3f] m/s^2\n', min(uz), max(uz));
    fprintf('vz range: [%.3f, %.3f] m/s^2\n', min(vz), max(vz));
end

%% Test B: 1D + uz >= 0 (thrust only upward)
fprintf('\n--- Test B: 1D + uz >= 0 ---\n');
cvx_begin quiet
    variable z2(1,N+1)
    variable vz2(1,N+1)
    variable uz2(1,N+1)
    minimize(0)
    subject to
        z2(1)==51; vz2(1)==0;
        z2(N+1)==1; vz2(N+1)==0;
        for k=1:N
            z2(k+1)==z2(k)+dt/2*(vz2(k)+vz2(k+1));
            vz2(k+1)==vz2(k)+dt/2*((uz2(k)-g0)+(uz2(k+1)-g0));
        end
        uz2>=0;
        z2>=0;
cvx_end
fprintf('Status: %s\n', cvx_status);
if strcmp(cvx_status,'Solved')||strcmp(cvx_status,'Inaccurate/Solved')
    fprintf('uz range: [%.3f, %.3f] m/s^2\n', min(uz2), max(uz2));
end

%% Test C: 1D + uz >= 0 + uz <= u_max
fprintf('\n--- Test C: 1D + 0 <= uz <= u_max ---\n');
cvx_begin quiet
    variable z3(1,N+1)
    variable vz3(1,N+1)
    variable uz3(1,N+1)
    minimize(0)
    subject to
        z3(1)==51; vz3(1)==0;
        z3(N+1)==1; vz3(N+1)==0;
        for k=1:N
            z3(k+1)==z3(k)+dt/2*(vz3(k)+vz3(k+1));
            vz3(k+1)==vz3(k)+dt/2*((uz3(k)-g0)+(uz3(k+1)-g0));
        end
        uz3>=0; uz3<=u_max;
        z3>=0;
cvx_end
fprintf('Status: %s\n', cvx_status);
if strcmp(cvx_status,'Solved')||strcmp(cvx_status,'Inaccurate/Solved')
    fprintf('uz range: [%.3f, %.3f] m/s^2\n', min(uz3), max(uz3));
    fprintf('vz range: [%.3f, %.3f] m/s^2\n', min(vz3), max(vz3));
end

%% Test D: 1D + vz rate limit
fprintf('\n--- Test D: 1D + 0<=uz<=u_max + |vz|<=7 ---\n');
cvx_begin quiet
    variable z4(1,N+1)
    variable vz4(1,N+1)
    variable uz4(1,N+1)
    minimize(0)
    subject to
        z4(1)==51; vz4(1)==0;
        z4(N+1)==1; vz4(N+1)==0;
        for k=1:N
            z4(k+1)==z4(k)+dt/2*(vz4(k)+vz4(k+1));
            vz4(k+1)==vz4(k)+dt/2*((uz4(k)-g0)+(uz4(k+1)-g0));
        end
        uz4>=0; uz4<=u_max;
        vz4>=-7; vz4<=7;
        z4>=0;
cvx_end
fprintf('Status: %s\n', cvx_status);
if strcmp(cvx_status,'Solved')||strcmp(cvx_status,'Inaccurate/Solved')
    fprintf('uz range: [%.3f, %.3f] m/s^2\n', min(uz4), max(uz4));
    fprintf('vz range: [%.3f, %.3f] m/s^2\n', min(vz4), max(vz4));
end

%% Test E: Same as D but tf=10 (shorter)
fprintf('\n--- Test E: same as D, tf=10s ---\n');
N2=40; tf2=10; dt2=tf2/N2;
cvx_begin quiet
    variable z5(1,N2+1)
    variable vz5(1,N2+1)
    variable uz5(1,N2+1)
    minimize(0)
    subject to
        z5(1)==51; vz5(1)==0;
        z5(N2+1)==1; vz5(N2+1)==0;
        for k=1:N2
            z5(k+1)==z5(k)+dt2/2*(vz5(k)+vz5(k+1));
            vz5(k+1)==vz5(k)+dt2/2*((uz5(k)-g0)+(uz5(k+1)-g0));
        end
        uz5>=0; uz5<=u_max;
        vz5>=-7; vz5<=7;
        z5>=0;
cvx_end
fprintf('Status: %s\n', cvx_status);
if strcmp(cvx_status,'Solved')||strcmp(cvx_status,'Inaccurate/Solved')
    fprintf('uz range: [%.3f, %.3f] m/s^2\n', min(uz5), max(uz5));
    fprintf('vz range: [%.3f, %.3f] m/s^2\n', min(vz5), max(vz5));
end