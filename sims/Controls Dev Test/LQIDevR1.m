 clear all; close all; clc;


% states = [x  3]

%xmaxs = [1 1 1 1 1 1 1e-3 1e-3 1e-3 1e-4 1e-4 1e-4];
%umaxs = [20 pi/16 pi/16 5];

%Qdiags = 1./xmaxs.^2;
%Rdiags = 1./umaxs.^2;


Qdiags = [0 0 0 10 10 10 1e2 1e2 1e2 1e6 1e6 1e6];
QIdiags =[10 10 500];
QI = diag(horzcat(Qdiags,QIdiags));
Q = diag(Qdiags);

% T delta rcs
Rdiags =[100,1e8,1e8,1e4];

R = diag(Rdiags);


params = struct();
params.g = 9.81;

params.d  = .2;
params.r_m = 0.125;
params.n_rcs = 1;
params.I = [1 0 0; 0 1 0; 0 0 1];

m = 1;

q = [sqrt(2)/2;0;sqrt(2)/2;0];
%q=[0 0 0 0]';
p = qToMRP(q)

x = [0;0;0;0;0;0;0;0;0;p];
u  = [10;0;0;0];

[A, B] = HopperLinearization_lqi(x, u, params, m);



%%

[Atilde,Btilde] = extendedDynamics(A,B);



CrossTerm = zeros(15,4);
%sysI = sys(Atilde,Btilde,CI,DI);

KlqI = lqr(Atilde,Btilde,QI,R,CrossTerm)
% KlqI(3,:) = KlqI(2,:);
% KlqI(3,2)= KlqI(2,1);
% KlqI(3,5)= KlqI(2,4);
% KlqI(3,9) = KlqI(2,8);
% KlqI(3,10) = KlqI(2,11);
% KlqI(3,14) = KlqI(2,13);
%KlqI(3,12) = -KlqI(2,11);
 %Klqi = lqi(sys,Q,R,CrossTerm);

%KlqI(:,10:12)