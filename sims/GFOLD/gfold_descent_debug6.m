%% Debug 6 — Bridge from working 1D to failing 3D
clear; clc;

g0=9.80665; m0=127.14;
Tmax=1.10*500*4.44822;
u_max=Tmax/m0;  % use m0 (not m_dry) — same as v2 script
r0=[0;0;51]; v0=[0;0;0]; rf=[0;0;1]; vf=[0;0;0];
g_vec=[0;0;-g0];
N=40; tf=20; dt=tf/N;

%% Test A: 3D kinematics, T vector, no bounds — mirror of debug5 Test A
fprintf('--- Test A: 3D free thrust vector ---\n');
cvx_begin quiet
    variable r(3,N+1); variable v(3,N+1); variable T(3,N+1)
    minimize(0)
    subject to
        r(:,1)==r0; v(:,1)==v0; r(:,N+1)==rf; v(:,N+1)==vf;
        for k=1:N
            r(:,k+1)==r(:,k)+dt/2*(v(:,k)+v(:,k+1));
            v(:,k+1)==v(:,k)+dt/2*((T(:,k)/m0+g_vec)+(T(:,k+1)/m0+g_vec));
        end
        r(3,:)>=0;
cvx_end
fprintf('Status: %s\n', cvx_status);
if strcmp(cvx_status,'Solved')||strcmp(cvx_status,'Inaccurate/Solved')
    fprintf('T(3) range: [%.2f, %.2f] N\n', min(T(3,:)), max(T(3,:)));
end

%% Test B: 3D + Tm magnitude slack + SOC only (no lower bound)
fprintf('\n--- Test B: 3D + SOC norm(T)<=Tm, Tm<=Tmax ---\n');
cvx_begin quiet
    variable r(3,N+1); variable v(3,N+1)
    variable T(3,N+1); variable Tm(1,N+1)
    minimize(0)
    subject to
        r(:,1)==r0; v(:,1)==v0; r(:,N+1)==rf; v(:,N+1)==vf;
        for k=1:N
            r(:,k+1)==r(:,k)+dt/2*(v(:,k)+v(:,k+1));
            v(:,k+1)==v(:,k)+dt/2*((T(:,k)/m0+g_vec)+(T(:,k+1)/m0+g_vec));
        end
        for k=1:N+1; norm(T(:,k),2)<=Tm(k); Tm(k)<=Tmax; end
        r(3,:)>=0; Tm>=0;
cvx_end
fprintf('Status: %s\n', cvx_status);

%% Test C: 3D + SOC + pointing T(3)>=Tm*cos(15deg)
fprintf('\n--- Test C: + pointing T(3)>=Tm*cos(15deg) ---\n');
cvx_begin quiet
    variable r(3,N+1); variable v(3,N+1)
    variable T(3,N+1); variable Tm(1,N+1)
    minimize(0)
    subject to
        r(:,1)==r0; v(:,1)==v0; r(:,N+1)==rf; v(:,N+1)==vf;
        for k=1:N
            r(:,k+1)==r(:,k)+dt/2*(v(:,k)+v(:,k+1));
            v(:,k+1)==v(:,k)+dt/2*((T(:,k)/m0+g_vec)+(T(:,k+1)/m0+g_vec));
        end
        for k=1:N+1
            norm(T(:,k),2)<=Tm(k); Tm(k)<=Tmax;
            T(3,k)>=Tm(k)*cos(deg2rad(15));
        end
        r(3,:)>=0; Tm>=0;
cvx_end
fprintf('Status: %s\n', cvx_status);

%% Test D: + vz limit
fprintf('\n--- Test D: + |vz|<=7 ---\n');
cvx_begin quiet
    variable r(3,N+1); variable v(3,N+1)
    variable T(3,N+1); variable Tm(1,N+1)
    minimize(0)
    subject to
        r(:,1)==r0; v(:,1)==v0; r(:,N+1)==rf; v(:,N+1)==vf;
        for k=1:N
            r(:,k+1)==r(:,k)+dt/2*(v(:,k)+v(:,k+1));
            v(:,k+1)==v(:,k)+dt/2*((T(:,k)/m0+g_vec)+(T(:,k+1)/m0+g_vec));
        end
        for k=1:N+1
            norm(T(:,k),2)<=Tm(k); Tm(k)<=Tmax;
            T(3,k)>=Tm(k)*cos(deg2rad(15));
            v(3,k)>=-7; v(3,k)<=7;
        end
        r(3,:)>=0; Tm>=0;
cvx_end
fprintf('Status: %s\n', cvx_status);

%% Test E: + mass tracking
fprintf('\n--- Test E: + mass tracking ---\n');
Isp=225;
cvx_begin quiet
    variable r(3,N+1); variable v(3,N+1)
    variable T(3,N+1); variable Tm(1,N+1); variable m(1,N+1)
    minimize(sum(Tm)*dt)
    subject to
        r(:,1)==r0; v(:,1)==v0; r(:,N+1)==rf; v(:,N+1)==vf; m(1)==m0;
        for k=1:N
            r(:,k+1)==r(:,k)+dt/2*(v(:,k)+v(:,k+1));
            v(:,k+1)==v(:,k)+dt/2*((T(:,k)/m0+g_vec)+(T(:,k+1)/m0+g_vec));
            m(k+1)==m(k)-dt/(2*Isp*g0)*(Tm(k)+Tm(k+1));
        end
        for k=1:N+1
            norm(T(:,k),2)<=Tm(k); Tm(k)<=Tmax;
            T(3,k)>=Tm(k)*cos(deg2rad(15));
            v(3,k)>=-7; v(3,k)<=7;
        end
        m>=119.90; m<=m0;
        r(3,:)>=0; Tm>=0;
cvx_end
fprintf('Status: %s\n', cvx_status);
if strcmp(cvx_status,'Solved')||strcmp(cvx_status,'Inaccurate/Solved')
    fprintf('Fuel used: %.2f kg\n', m0-m(end));
    fprintf('Max |vz|: %.3f m/s\n', max(abs(v(3,:))));
    fprintf('Final pos: [%.3f %.3f %.3f] m\n', r(1,end),r(2,end),r(3,end));
end

%% Test F: Drop m<=m0, remove redundant upper bound
fprintf('\n--- Test F: mass tracking, drop m<=m0 ---\n');
Isp=225;
cvx_begin quiet
    variable r(3,N+1); variable v(3,N+1)
    variable T(3,N+1); variable Tm(1,N+1); variable m(1,N+1)
    minimize(sum(Tm)*dt)
    subject to
        r(:,1)==r0; v(:,1)==v0; r(:,N+1)==rf; v(:,N+1)==vf; m(1)==m0;
        for k=1:N
            r(:,k+1)==r(:,k)+dt/2*(v(:,k)+v(:,k+1));
            v(:,k+1)==v(:,k)+dt/2*((T(:,k)/m0+g_vec)+(T(:,k+1)/m0+g_vec));
            m(k+1)==m(k)-dt/(2*Isp*g0)*(Tm(k)+Tm(k+1));
        end
        for k=1:N+1
            norm(T(:,k),2)<=Tm(k); Tm(k)<=Tmax;
            T(3,k)>=Tm(k)*cos(deg2rad(15));
            v(3,k)>=-7; v(3,k)<=7;
        end
        m>=119.90;          % only lower bound — m<=m0 is redundant
        r(3,:)>=0; Tm>=0;
cvx_end
fprintf('Status: %s\n', cvx_status);
if strcmp(cvx_status,'Solved')||strcmp(cvx_status,'Inaccurate/Solved')
    fprintf('Fuel used: %.2f kg\n', m0-m(end));
    fprintf('Max |vz|: %.3f m/s\n', max(abs(v(3,:))));
    fprintf('Final pos: [%.3f %.3f %.3f] m\n', r(1,end),r(2,end),r(3,end));
    fprintf('Tm range: [%.1f, %.1f] N\n', min(Tm), max(Tm));
end

%% Test G: Same but cap Tmax by prop budget
fprintf('\n--- Test G: Tmax capped by prop budget ---\n');
% max avg thrust such that total impulse <= prop_remain * Isp * g0
prop_remain=7.24;
Tmax_budget = prop_remain*Isp*g0/tf;  % average thrust cap
fprintf('Tmax_budget = %.1f N\n', Tmax_budget);
cvx_begin quiet
    variable r(3,N+1); variable v(3,N+1)
    variable T(3,N+1); variable Tm(1,N+1); variable m(1,N+1)
    minimize(sum(Tm)*dt)
    subject to
        r(:,1)==r0; v(:,1)==v0; r(:,N+1)==rf; v(:,N+1)==vf; m(1)==m0;
        for k=1:N
            r(:,k+1)==r(:,k)+dt/2*(v(:,k)+v(:,k+1));
            v(:,k+1)==v(:,k)+dt/2*((T(:,k)/m0+g_vec)+(T(:,k+1)/m0+g_vec));
            m(k+1)==m(k)-dt/(2*Isp*g0)*(Tm(k)+Tm(k+1));
        end
        for k=1:N+1
            norm(T(:,k),2)<=Tm(k);
            Tm(k)<=Tmax;           % hardware limit
            T(3,k)>=Tm(k)*cos(deg2rad(15));
            v(3,k)>=-7; v(3,k)<=7;
        end
        m>=119.90;
        % Prop budget as global constraint
        sum(Tm)*dt <= prop_remain*Isp*g0;
        r(3,:)>=0; Tm>=0;
cvx_end
fprintf('Status: %s\n', cvx_status);
if strcmp(cvx_status,'Solved')||strcmp(cvx_status,'Inaccurate/Solved')
    fprintf('Fuel used: %.2f kg\n', m0-m(end));
    fprintf('Max |vz|: %.3f m/s\n', max(abs(v(3,:))));
    fprintf('Final pos: [%.3f %.3f %.3f] m\n', r(1,end),r(2,end),r(3,end));
end

%% Test H: No mass variable at all — post-process mass from Tm solution
fprintf('\n--- Test H: no mass variable, minimize sum(Tm), post-process mass ---\n');
Isp_val=225;
cvx_begin quiet
    variable r(3,N+1); variable v(3,N+1)
    variable T(3,N+1); variable Tm(1,N+1)
    minimize(sum(Tm)*dt)
    subject to
        r(:,1)==r0; v(:,1)==v0; r(:,N+1)==rf; v(:,N+1)==vf;
        for k=1:N
            r(:,k+1)==r(:,k)+dt/2*(v(:,k)+v(:,k+1));
            v(:,k+1)==v(:,k)+dt/2*((T(:,k)/m0+g_vec)+(T(:,k+1)/m0+g_vec));
        end
        for k=1:N+1
            norm(T(:,k),2)<=Tm(k);
            Tm(k)<=Tmax;
            T(3,k)>=Tm(k)*cos(deg2rad(15));
            v(3,k)>=-7; v(3,k)<=7;
        end
        % Prop budget as single linear constraint (no mass chain)
        sum(Tm)*dt <= 7.24*Isp_val*g0;
        r(3,:)>=0; Tm>=0;
cvx_end
fprintf('Status: %s\n', cvx_status);
if strcmp(cvx_status,'Solved')||strcmp(cvx_status,'Inaccurate/Solved')
    % Post-process mass
    m_post = zeros(1,N+1); m_post(1)=m0;
    for k=1:N
        m_post(k+1)=m_post(k)-dt/(2*Isp_val*g0)*(Tm(k)+Tm(k+1));
    end
    fprintf('Fuel used: %.2f kg\n', m0-m_post(end));
    fprintf('Min mass: %.2f kg (dry=%.2f)\n', min(m_post), 119.90);
    fprintf('Max |vz|: %.3f m/s\n', max(abs(v(3,:))));
    fprintf('Final pos: [%.3f %.3f %.3f] m\n', r(1,end),r(2,end),r(3,end));
    fprintf('Tm range: [%.1f %.1f] N\n', min(Tm), max(Tm));
end

%% Test I: Same as H but with maximize(-sum(Tm)) instead of minimize
fprintf('\n--- Test I: same as H, feasibility only (minimize 0) ---\n');
cvx_begin quiet
    variable r(3,N+1); variable v(3,N+1)
    variable T(3,N+1); variable Tm(1,N+1)
    minimize(0)
    subject to
        r(:,1)==r0; v(:,1)==v0; r(:,N+1)==rf; v(:,N+1)==vf;
        for k=1:N
            r(:,k+1)==r(:,k)+dt/2*(v(:,k)+v(:,k+1));
            v(:,k+1)==v(:,k)+dt/2*((T(:,k)/m0+g_vec)+(T(:,k+1)/m0+g_vec));
        end
        for k=1:N+1
            norm(T(:,k),2)<=Tm(k);
            Tm(k)<=Tmax;
            T(3,k)>=Tm(k)*cos(deg2rad(15));
            v(3,k)>=-7; v(3,k)<=7;
        end
        sum(Tm)*dt <= 7.24*Isp_val*g0;
        r(3,:)>=0; Tm>=0;
cvx_end
fprintf('Status: %s\n', cvx_status);

%% Test J: Remove prop budget constraint too
fprintf('\n--- Test J: no mass, no prop budget, just SOC+pointing+vz+BCs ---\n');
cvx_begin quiet
    variable r(3,N+1); variable v(3,N+1)
    variable T(3,N+1); variable Tm(1,N+1)
    minimize(0)
    subject to
        r(:,1)==r0; v(:,1)==v0; r(:,N+1)==rf; v(:,N+1)==vf;
        for k=1:N
            r(:,k+1)==r(:,k)+dt/2*(v(:,k)+v(:,k+1));
            v(:,k+1)==v(:,k)+dt/2*((T(:,k)/m0+g_vec)+(T(:,k+1)/m0+g_vec));
        end
        for k=1:N+1
            norm(T(:,k),2)<=Tm(k);
            Tm(k)<=Tmax;
            T(3,k)>=Tm(k)*cos(deg2rad(15));
            v(3,k)>=-7; v(3,k)<=7;
        end
        r(3,:)>=0; Tm>=0;
cvx_end
fprintf('Status: %s\n', cvx_status);
if strcmp(cvx_status,'Solved')||strcmp(cvx_status,'Inaccurate/Solved')
    fprintf('This matches Test D — confirming baseline\n');
    fprintf('Tm range: [%.1f %.1f] N\n', min(Tm), max(Tm));
    fprintf('sum(Tm)*dt = %.1f N*s  (prop equiv = %.3f kg)\n', ...
        sum(Tm)*dt, sum(Tm)*dt/(Isp_val*g0));
end

%% Test K: Shorter tf — force ballistic + braking profile
% Free-fall from 51->1m takes ~3.2s minimum
% With vz=7 cap, need to start braking before 7 m/s
% Try tf=8s — enough for fall + braking, not so long solver idles
fprintf('\n--- Test K: tf=8s, no prop budget ---\n');
Isp_val=225; N3=40; tf3=8; dt3=tf3/N3;
cvx_begin quiet
    variable r(3,N3+1); variable v(3,N3+1)
    variable T(3,N3+1); variable Tm(1,N3+1)
    minimize(sum(Tm)*dt3)
    subject to
        r(:,1)==r0; v(:,1)==v0; r(:,N3+1)==rf; v(:,N3+1)==vf;
        for k=1:N3
            r(:,k+1)==r(:,k)+dt3/2*(v(:,k)+v(:,k+1));
            v(:,k+1)==v(:,k)+dt3/2*((T(:,k)/m0+g_vec)+(T(:,k+1)/m0+g_vec));
        end
        for k=1:N3+1
            norm(T(:,k),2)<=Tm(k); Tm(k)<=Tmax;
            T(3,k)>=Tm(k)*cos(deg2rad(15));
            v(3,k)>=-7; v(3,k)<=7;
        end
        r(3,:)>=0; Tm>=0;
cvx_end
fprintf('Status: %s\n', cvx_status);
if strcmp(cvx_status,'Solved')||strcmp(cvx_status,'Inaccurate/Solved')
    fprintf('sum(Tm)*dt = %.1f N*s (prop equiv = %.3f kg)\n', ...
        sum(Tm)*dt3, sum(Tm)*dt3/(Isp_val*g0));
    fprintf('Tm range: [%.1f %.1f] N\n', min(Tm), max(Tm));
    fprintf('vz range: [%.3f %.3f] m/s\n', min(v(3,:)), max(v(3,:)));
end

%% Test L: tf=8s + prop budget constraint
fprintf('\n--- Test L: tf=8s + prop budget 7.24kg ---\n');
cvx_begin quiet
    variable r(3,N3+1); variable v(3,N3+1)
    variable T(3,N3+1); variable Tm(1,N3+1)
    minimize(sum(Tm)*dt3)
    subject to
        r(:,1)==r0; v(:,1)==v0; r(:,N3+1)==rf; v(:,N3+1)==vf;
        for k=1:N3
            r(:,k+1)==r(:,k)+dt3/2*(v(:,k)+v(:,k+1));
            v(:,k+1)==v(:,k)+dt3/2*((T(:,k)/m0+g_vec)+(T(:,k+1)/m0+g_vec));
        end
        for k=1:N3+1
            norm(T(:,k),2)<=Tm(k); Tm(k)<=Tmax;
            T(3,k)>=Tm(k)*cos(deg2rad(15));
            v(3,k)>=-7; v(3,k)<=7;
        end
        sum(Tm)*dt3 <= 7.24*Isp_val*g0;
        r(3,:)>=0; Tm>=0;
cvx_end
fprintf('Status: %s\n', cvx_status);
if strcmp(cvx_status,'Solved')||strcmp(cvx_status,'Inaccurate/Solved')
    fprintf('Fuel used: %.3f kg\n', sum(Tm)*dt3/(Isp_val*g0));
    fprintf('vz range: [%.3f %.3f] m/s\n', min(v(3,:)), max(v(3,:)));
    fprintf('Final pos: [%.3f %.3f %.3f] m\n', r(1,end),r(2,end),r(3,end));
end