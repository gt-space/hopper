clear all; close all; clc;

global g;
global b;
global m; 
global h;
global a;

% Define constants and parameters
g = 9.81; % acceleration due to gravity [m/s^2]
b = 0.5;    % thrust weighting
m = 10;   % mass [kg]
h = -53;  % target height [m]
a=200;

X0_ascent =[3.4; 19.6; 5.7];

k=3;
tf0 = sqrt(-2*h/g); 
a10 = 0;
a20 = (k*g)/(1+k);      
X0_descent  = [a10; a20; tf0];

%X0_descent = [0; 0; sqrt(-2*h/g) ];

X = fsolve(@OCP,X0_ascent);
X2 = fsolve(@OCP_Descent,X0_descent);

a1_ascent = X(1);
a2_ascent = X(2);
tf_ascent = X(3);
t_ascent = linspace(0,tf_ascent,100);

z_ascent = 0.5*g.*t_ascent.^2 + (1/6)*a1_ascent.*t_ascent.^3 - 0.5*a2_ascent.*t_ascent.^2;
v_ascent = g.*t_ascent + 0.5*a1_ascent.*t_ascent.^2 - a2_ascent.*t_ascent;

lv =  m^2*(-a1_ascent.*t_ascent+a2_ascent)*b;
u_ascent = lv./(m);

%% Descent

a1_descent = X2(1);
a2_descent = X2(2);
tf_descent = X2(3);

t_descent = linspace(0,tf_descent,100);

z_descent = 0.5*g.*t_descent.^2 + (1/6)*a1_descent.*t_descent.^3 - 0.5*a2_descent.*t_descent.^2 + h;

% FIX (explicitly wrong): t_ascent -> t_descent
v_descent = g.*t_descent + 0.5*a1_descent.*t_descent.^2 - a2_descent.*t_descent;

lv_descent =  m^2*(-a1_descent.*t_descent+a2_descent)*b;
u_descent = lv_descent./(m);

%% Concatenate ascent + descent, cutting the repeated junction sample EARLY
% Drop the first descent point so time samples are unique at the join.
t  = horzcat(t_ascent, (t_descent(2:end) + tf_ascent));
z  = horzcat(z_ascent, z_descent(2:end));
vz = horzcat(v_ascent, v_descent(2:end));
T_profile = horzcat(u_ascent, u_descent(2:end));

figure
plot(t,z);
xlabel('t [s]');
ylabel('z [m]');
grid on
title('z trajectory (+z down)');

figure
plot(t,T_profile)
xlabel('t [s]');
ylabel('T [N]');
grid on
title('Thrust Profile');

%%

stateMaxs = [100; 100;100;0.5; pi; 4];
controlMaxs = [100; pi/8];

tf = t(end);
t0=0;

% 1/A(ii) = (tf-t0)*max acceptable value of x_i(t)^2
Q = [  (1/stateMaxs(1))^2 0 0 0 0 0; 
       0 (1/stateMaxs(2))^2 0 0 0 0; 
       0 0 (1/stateMaxs(3))^2 0 0 0; 
       0 0 0 (1/stateMaxs(4))^2 0 0; 
       0 0 0 0 (1/stateMaxs(5))^2 0;
       0 0 0 0 0 (1/stateMaxs(6))^2].*(tf-t0);

%  1/B(ii) = (tf-t0)*max acceptable value of u_i(t)^2
R = [(1/controlMaxs(1)^2) 0
     0 1/controlMaxs(2)^2].*(tf-t0);

% x z vx vx theta omega
finalStateMaxs = [10; -80; 1; 1; pi; .1];

theta = ones(1,length(t)).*pi/2;
delta = zeros(1,length(t));

% Linearize along the (now unique) full trajectory
for i = 1:length(t)
    [A(:,:,i),B(:,:,i)]  = HopperLinearization_3DOF(theta(i),T_profile(i),delta(i));
end

%% target trajectory (keep same length as t)
x_trajectory = zeros(1,length(t));
z_trajectory = z;
vx_trajectory = zeros(1,length(t));
vz_trajectory = vz;
theta_trajectory = (pi/2).*ones(1,length(t));
q_trajectory = zeros(1,length(t));

Target_Trajectory = [x_trajectory;z_trajectory;vx_trajectory;vz_trajectory;theta_trajectory;q_trajectory]';

% Use the cleaned mission time grid directly
tgrid = t;

% Call gains using consistent (unique) grid + consistent trajectory
GainTable = TrajectoryFollowingGains(finalStateMaxs,stateMaxs,controlMaxs,Target_Trajectory,A,B,tgrid);

%%

K1 = GainTable.K1;
K2grid = GainTable.K2.';   % keep your Simulink convention (N×2)

unom = [T_profile(:), zeros(length(T_profile),1)];  % N×2

% Flatten K1 for Simulink (N×12)
K1flat = zeros(length(tgrid),12);
for i = 1:length(tgrid)
    K1flat(i,:) = reshape(K1(:,:,i).',1,12);
end

%% --- your OCP functions unchanged below ---

function F = OCP(X)
global g; global b; global m; global h;

a1 = X(1);
a2 = X(2);
tf = X(3);

F = zeros(3,1);

% vf = 0
F(1) = (g*tf)+(0.5*(a1*tf^2))-(a2*tf);

% zf = h
F(2) = (0.5*g*tf^2)+((1/6)*(a1*tf^3))-(0.5*a2*tf^2)-h;

% dPhi/dtf +Htf = 0
nu2 = m^2*(-a1*tf+a2)/b;
F(3) = nu2*g - 0.5*nu2^2/(b*m^2);
end

function F = OCP_Descent(X)
global g; global b; global m; global h; global a;

a1 = X(1);
a2 = X(2);
tf = X(3);

F = zeros(3,1);

% zf = 0
F(1) = (0.5*g*tf^2)+((1/6)*(a1*tf^3))-(0.5*a2*tf^2)+h;

vf = (g*tf)+(0.5*(a1*tf^2))-(a2*tf);
lvf = (b*m^2)*(-a1*tf+a2);

F(2) = lvf-a*vf;

% dPhi/dtf +Htf = 0
uf = lvf/(b*m);
lzf = a1*b*m^2;
F(3) = (0.5*b*uf^2)+(lzf*vf)+lvf*(g-uf/m);
end
