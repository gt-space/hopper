clear
clc
close all

%% GENERATE TRAJECTORY

[Trajectory,unom,tgrid,m_profile,I_profile] = Trajectory3();

N = length(tgrid);

nx = 13;
nu = 4;

%% COST MATRICES

Sf = diag([1 1 1e10 1 1 1e6 1 1 1 1000 1000 1000 1000]);
Q  = diag([1 1 1e10 1 1 10 1 1 1 100 100 100 100]);
R  = diag([100 1000 1000 1000]);

%% PARAMETERS

params = struct();

params.g = 9.81;
params.d = 1.0;
params.r_m = 0.125;

%% PREALLOCATE

A=zeros(nx,nx,N);
B=zeros(nx,nu,N);
Fref=zeros(N,nx);

%% LINEARIZATION LOOP

for i=1:N

    params.m = m_profile(i);
    params.I = I_profile(:,:,i);

    x = Trajectory(i,:)';
    u = unom(i,:)';

    [Ai,Bi,f0] = HopperLinearization_6DOF_inertial(x,u,params);
    

    A(:,:,i)=Ai;
    B(:,:,i)=Bi;
    Fref(i,:)=f0';

end

%% COMPUTE LQR

GainTable = TrajectoryFollowingGains2( ...
Sf,Q,R,Trajectory,A,B,tgrid,unom,params,Fref);

K1 = GainTable.K1;
K2grid = GainTable.K2';

%% FLATTEN K1 FOR SIMULINK

K1flat=zeros(N,52);

for i=1:N
    K1flat(i,:) = reshape(K1(:,:,i).',1,52);
end

disp("LQR Gains Generated")