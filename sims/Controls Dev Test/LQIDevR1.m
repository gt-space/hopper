 clear all; close all; clc;


% states = [x y z vx vy vz P Q R p1 p2 p3]


Qdiags = [1e-4 1e-4 1e-4 1 1 1 1 1 1 1 1e-9 1];
QIdiags =[500 500 500 500];
QI = diag(horzcat(Qdiags,QIdiags));
Q = diag(Qdiags);

% T delta rcs
Rdiags =[1,1000,1000,100];

R = diag(Rdiags);


params = struct();
params.g = 9.80665;

params.d  = .2;
params.r_m = 0.125;
params.n_rcs = 1;
params.I = eye(3);

m = 1;

q = [sqrt(2)/2;0;sqrt(2)/2;0];
p = qToMRP(q)

x = [0;0;0;0;0;0;0;0;0;p];
u  = [10;0;0;0];

[A, B] = HopperLinearization_lqi(x, u, params, m);



%GainTable = TrajectoryFollowingGains3(Sf,Q,R,Target_Trajectory,A,B,tgrid,unom,params,Fref);

C = eye(12);
D = zeros(12,4);
sys = ss(A,B,C,D);



Klqr = lqr(sys,Q,R)

%%

[Atilde,Btilde] = extendedDynamics(A,B);
CI =eye(16);
DI = zeros(16,4);
CrossTerm = zeros(16,4);
%sysI = sys(Atilde,Btilde,CI,DI);

KlqI = lqr(Atilde,Btilde,QI,R,CrossTerm)

 %Klqi = lqi(sys,Q,R,CrossTerm);

