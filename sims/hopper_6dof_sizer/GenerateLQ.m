
%% ================= GENERATE TRAJECTORY =================

[Trajectory, unom, tgrid] = Trajectory3();

N  = length(tgrid);
nx = 13;
nu = 4;

%% ================= COST MATRICES =================

Sf = diag([ ...
1 ...
1 ...
1e10 ...
1 ...
1 ...
1e6 ...
1 ...
1 ...
1 ...
1000 ...
1000 ...
1000 ...
1000 ]);

Q = diag([ ...
1 ...
1 ...
1e8 ...
1 ...
1 ...
10 ...
1 ...
1 ...
1 ...
100 ...
100 ...
100 ...
100 ]);

R = diag([ ...
2e6 ...
1000 ...
1000 ...
1000 ]);

%% ================= VEHICLE PARAMETERS =================

params = struct();

params.g = 9.80665;
params.m = 114;

params.d   = 1.0;
params.r_m = 0.125;

Ix = 16.5;
Iy = 16.5;
Iz = 16.5;

params.I = diag([Ix Iy Iz]);

%% ================= LINEARIZE ALONG TRAJECTORY =================

A = zeros(nx,nx,N);
B = zeros(nx,nu,N);
Fref = zeros(N,nx);

for i = 1:N

    x = Trajectory(i,:)';
    u = unom(i,:)';

    [Ai,Bi,f0] = HopperLinearization_6DOF(x,u,params);

    A(:,:,i) = Ai;
    B(:,:,i) = Bi;
    Fref(i,:) = f0';

end

%% ================= COMPUTE TVLQR GAINS =================

GainTable = TrajectoryFollowingGains2( ...
    Sf,Q,R,Trajectory,A,B,tgrid,unom,params,Fref);

K1 = GainTable.K1;
K2 = GainTable.K2';

%% ================= FLATTEN FOR SIMULINK =================

K1flat = zeros(N,52);

for i = 1:N
    K1flat(i,:) = reshape(K1(:,:,i).',1,52);
end

%% ================= EXPORT TIMESERIES =================

z_ref_ts  = timeseries(Trajectory(:,3),tgrid);
vz_ref_ts = timeseries(Trajectory(:,6),tgrid);
T_ref_ts  = timeseries(unom(:,1),tgrid);

disp("TVLQR Gain Generation Complete")