%% Debug 4 — Test with zero Tmin (free coast) and relaxed pointing
clear; clc;

m0=127.14; m_dry=119.90; g0=9.80665; Isp=225; alpha=1/(Isp*g0);
T_nom=500*4.44822; Tmax=1.10*T_nom;
u_max=Tmax/m_dry;
r0=[0;0;51]; v0=[0;0;0]; rf=[0;0;1]; vf=[0;0;0];
vz_max=7; g_vec=[0;0;-g0];
N=40; tf=20; dt=tf/N; t=linspace(0,tf,N+1);

%% Test A: No Tmin (free coast), no pointing, full z chain
fprintf('--- Test A: Tmin=0, no pointing ---\n');
cvx_begin quiet
    variable rA(3,N+1); variable vA(3,N+1)
    variable uA(3,N+1); variable sA(1,N+1); variable zA(1,N+1)
    maximize(zA(N+1))
    subject to
        rA(:,1)==r0; vA(:,1)==v0; rA(:,N+1)==rf; vA(:,N+1)==vf; zA(1)==log(m0);
        for k=1:N
            rA(:,k+1)==rA(:,k)+dt/2*(vA(:,k)+vA(:,k+1));
            vA(:,k+1)==vA(:,k)+dt/2*((uA(:,k)+g_vec)+(uA(:,k+1)+g_vec));
            zA(k+1)==zA(k)-dt/2*alpha*(sA(k)+sA(k+1));
        end
        for k=1:N+1
            norm(uA(:,k),2)<=sA(k);
            sA(k)<=u_max;
            vA(3,k)>=-vz_max; vA(3,k)<=vz_max;
        end
        zA>=log(m_dry); rA(3,:)>=0; sA>=0;
cvx_end
fprintf('Status: %s\n', cvx_status);
if strcmp(cvx_status,'Solved')||strcmp(cvx_status,'Inaccurate/Solved')
    uA_mag=sqrt(sum(uA.^2,1));
    angles_A=acosd(max(uA(3,:),0)./max(uA_mag,1e-8));
    fprintf('Fuel used: %.2f kg\n', m0-exp(zA(end))*m0/m0);
    fprintf('Thrust angle range: [%.1f, %.1f] deg\n', min(angles_A), max(angles_A));
    fprintf('Nodes with zero thrust (sA<0.01): %d\n', sum(sA<0.01));
end

%% Test B: Tmin=0, pointing 80deg (very loose)
fprintf('\n--- Test B: Tmin=0, theta_max=80deg ---\n');
cvx_begin quiet
    variable rB(3,N+1); variable vB(3,N+1)
    variable uB(3,N+1); variable sB(1,N+1); variable zB(1,N+1)
    maximize(zB(N+1))
    subject to
        rB(:,1)==r0; vB(:,1)==v0; rB(:,N+1)==rf; vB(:,N+1)==vf; zB(1)==log(m0);
        for k=1:N
            rB(:,k+1)==rB(:,k)+dt/2*(vB(:,k)+vB(:,k+1));
            vB(:,k+1)==vB(:,k)+dt/2*((uB(:,k)+g_vec)+(uB(:,k+1)+g_vec));
            zB(k+1)==zB(k)-dt/2*alpha*(sB(k)+sB(k+1));
        end
        for k=1:N+1
            norm(uB(:,k),2)<=sB(k);
            sB(k)<=u_max;
            uB(3,k)>=sB(k)*cos(deg2rad(80));
            vB(3,k)>=-vz_max; vB(3,k)<=vz_max;
        end
        zB>=log(m_dry); rB(3,:)>=0; sB>=0;
cvx_end
fprintf('Status: %s\n', cvx_status);
if strcmp(cvx_status,'Solved')||strcmp(cvx_status,'Inaccurate/Solved')
    fprintf('Fuel used: %.2f kg\n', m0-exp(zB(end)));
end

%% Test C: Tmin=0, pointing 45deg
fprintf('\n--- Test C: Tmin=0, theta_max=45deg ---\n');
cvx_begin quiet
    variable rC(3,N+1); variable vC(3,N+1)
    variable uC(3,N+1); variable sC(1,N+1); variable zC(1,N+1)
    maximize(zC(N+1))
    subject to
        rC(:,1)==r0; vC(:,1)==v0; rC(:,N+1)==rf; vC(:,N+1)==vf; zC(1)==log(m0);
        for k=1:N
            rC(:,k+1)==rC(:,k)+dt/2*(vC(:,k)+vC(:,k+1));
            vC(:,k+1)==vC(:,k)+dt/2*((uC(:,k)+g_vec)+(uC(:,k+1)+g_vec));
            zC(k+1)==zC(k)-dt/2*alpha*(sC(k)+sC(k+1));
        end
        for k=1:N+1
            norm(uC(:,k),2)<=sC(k);
            sC(k)<=u_max;
            uC(3,k)>=sC(k)*cos(deg2rad(45));
            vC(3,k)>=-vz_max; vC(3,k)<=vz_max;
        end
        zC>=log(m_dry); rC(3,:)>=0; sC>=0;
cvx_end
fprintf('Status: %s\n', cvx_status);
if strcmp(cvx_status,'Solved')||strcmp(cvx_status,'Inaccurate/Solved')
    fprintf('Fuel used: %.2f kg\n', m0-exp(zC(end)));
end

%% Test D: Tmin=0, pointing 10deg (original)
fprintf('\n--- Test D: Tmin=0, theta_max=10deg ---\n');
cvx_begin quiet
    variable rD(3,N+1); variable vD(3,N+1)
    variable uD(3,N+1); variable sD(1,N+1); variable zD(1,N+1)
    maximize(zD(N+1))
    subject to
        rD(:,1)==r0; vD(:,1)==v0; rD(:,N+1)==rf; vD(:,N+1)==vf; zD(1)==log(m0);
        for k=1:N
            rD(:,k+1)==rD(:,k)+dt/2*(vD(:,k)+vD(:,k+1));
            vD(:,k+1)==vD(:,k)+dt/2*((uD(:,k)+g_vec)+(uD(:,k+1)+g_vec));
            zD(k+1)==zD(k)-dt/2*alpha*(sD(k)+sD(k+1));
        end
        for k=1:N+1
            norm(uD(:,k),2)<=sD(k);
            sD(k)<=u_max;
            uD(3,k)>=sD(k)*cos(deg2rad(10));
            vD(3,k)>=-vz_max; vD(3,k)<=vz_max;
        end
        zD>=log(m_dry); rD(3,:)>=0; sD>=0;
cvx_end
fprintf('Status: %s\n', cvx_status);
if strcmp(cvx_status,'Solved')||strcmp(cvx_status,'Inaccurate/Solved')
    fprintf('Fuel used: %.2f kg\n', m0-exp(zD(end)));
end