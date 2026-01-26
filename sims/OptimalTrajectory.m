clear all; close all; clc;

global g;
global b;
global m; 
global h;
global a;

% Define constants and parameters
g = 9.81; % acceleration due to gravity [m/s^2]
b = 1;    % thrust weighting
m = 10;   % mass [kg]
h = -53;  % target height [m]
a=300;

X0_ascent =[3.4; 19.6; 5.7];
X0_descent = [0; 0; sqrt(-2*h/g) ];

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

z_descent = 0.5*g.*t_descent.^2 + (1/6)*a1_descent.*t_descent.^3 - 0.5*a2_descent.*t_descent.^2+h;
v_descent = g.*t_descent + 0.5*a1_descent.*t_descent.^2 - a2_descent.*t_ascent;

lv_descent =  m^2*(-a1_descent.*t_descent+a2_descent)*b;
u_descent = lv_descent./(m);

%% Plotting

t = horzcat(t_ascent,t_descent+tf_ascent);
z = horzcat(z_ascent,z_descent);
u = horzcat(u_ascent,u_descent);

figure
plot(t,z);
xlabel('t [s]');
ylabel('z [m]');
grid on
title('z trajectory (+z down)');


figure
plot(t,u)
xlabel('t [s]');
ylabel('T [N]');
grid on
title('Thrust Profile');


function F = OCP(X)

global g;
global b;
global m; 
global h;


a1 = X(1);
a2 = X(2);
tf = X(3);

F = zeros(3,1);

% I am doing this substitution because having large mass makes fsolve not converge due to 1/m^2 on c1
% also large changes in b cause similar issues
%a1 = c1/(b*m^2)
%a2 = c2/(b*m^2)

% vf = 0
F(1) = (g*tf)+(0.5*(a1*tf^2))-(a2*tf);

% zf = h
F(2) = (0.5*g*tf^2)+((1/6)*(a1*tf^3))-(0.5*a2*tf^2)-h;

% dPhi/dtf +Htf = 0
nu2 = m^2*(-a1*tf+a2)/b;
F(3) = nu2*g - 0.5*nu2^2/(b*m^2);
end

function F = OCP_Descent(X)

global g;
global b;
global m; 
global h;
global a;
a1 = X(1);
a2 = X(2);
tf = X(3);

F = zeros(3,1);

% I am doing this substitution because having large mass makes fsolve not converge due to 1/m^2 on c1
% also large changes in b cause similar issues
%a1 = c1/(b*m^2)
%a2 = c2/(b*m^2)

% zf = 0
F(1) = (0.5*g*tf^2)+((1/6)*(a1*tf^3))-(0.5*a2*tf^2)+h;

vf = (g*tf)+(0.5*(a1*tf^2))-(a2*tf);
lvf = (b*m^2)*(-a1*tf+a2);

F(2) =lvf-a*vf;

% dPhi/dtf +Htf = 0
uf = lvf/(b*m);
lzf = a1*b*m^2;
F(3) = (0.5*b*uf^2)+(lzf*vf)+lvf*(g-uf/m);
end