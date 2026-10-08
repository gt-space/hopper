 clear all; close all; clc;


% states = [x y z vx vy vz P Q R p1 p2 p3]

%xmaxs = [1 1 1 1 1 1 1e-3 1e-3 1e-3 1e-4 1e-4 1e-4];
%umaxs = [20 pi/16 pi/16 5];

%Qdiags = 1./xmaxs.^2;
%Rdiags = 1./umaxs.^2;


Qdiags = [0 0 0 30 30 30 1 1 1 10 10 10];
QIdiags =[10 10 500];
QI = diag(horzcat(Qdiags,QIdiags));
Q = diag(Qdiags);

% T delta rcs
Rdiags =[100,1e9,1e9,1];

R = diag(Rdiags);


params = struct();
params.g = 9.81;

params.d  = .2;
params.r_m = 0.125;
params.n_rcs = 1;
params.I = [1 0 0; 0 4 0; 0 0 4];

m = 1;

q = [sqrt(2)/2;0;sqrt(2)/2;0];
%q=[1 0 0 0]';
p = qToMRP(q)

x = [0;0;0;0;0;0;0;0;0;p];
u  = [10;1e-4;1e-4;0];

[A, B] = HopperLinearization_lqi(x, u, params, m);



%%

[Atilde,Btilde] = extendedDynamics(A,B);
CrossTerm = zeros(15,4);

%sysI = sys(Atilde,Btilde,CI,DI);

KlqI = lqr(Atilde,Btilde,QI,R,CrossTerm)
%Klqi = lqi(sys,Q,R,CrossTerm);
% KlqI(:,10:12)

CL = Atilde - Btilde*KlqI;
poles = eig(CL);
[V,D] = eig(CL);

poles
plot(real(poles), imag(poles),'rx')

tolerance = 1e-5; %Tune this

names = {'x','y','z','vx','vy','vz','P','Q','R','p1','p2','p3','ex','ey','ez'};
for k = 1:15
  vnormalized = V(:,k);
  idx    = find(vnormalized > tolerance);
  [vals, order] = sort(vnormalized(idx), 'descend');
  idx    = idx(order);
  %parts = arrayfun(@(j) sprintf('%s = %.3f', names{j}, v(j)), idx, 'UniformOutput', false);
  phi = rad2deg(angle(V(:,k)));
  parts = compose('%s(%+.0f°)', string(names(idx)).', phi(idx));
  fprintf('pole %8.3f%+8.3fi, phase= : %s\n', real(D(k,k)), imag(D(k,k)), strjoin(parts, ', '))
end