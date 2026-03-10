function [Trajectory,unom,tgrid,m_profile,I_profile] = Trajectory3()

close all

%% VEHICLE PARAMETERS

m0 = 114;
mf = 94;
g  = 9.81;

Tmin = 890;
Tmax = 2400;

h_target = 50;

t_ascend  = 10;
t_hover   = 5;
t_descend = 10;

dt = 0.02;

%% TIME GRID

t1 = 0:dt:t_ascend;
t2 = t1(end)+dt:dt:t1(end)+t_hover;
t3 = t2(end)+dt:dt:t2(end)+t_descend;

tgrid = [t1 t2 t3]';
N = length(tgrid);

tf = tgrid(end);

%% MASS PROFILE

m_profile = m0 - (m0-mf)*(tgrid/tf);

%% INERTIA PROFILE

Ix0 = 1;
Iy0 = 40;
Iz0 = 40;

I_profile = zeros(3,3,N);

for i=1:N

    scale = m_profile(i)/m0;

    I_profile(:,:,i) = diag([Ix0 Iy0 Iz0])*scale;

end

%% MINIMUM JERK ALTITUDE

z  = zeros(N,1);
vz = zeros(N,1);
az = zeros(N,1);

for i = 1:length(t1)

    s = t1(i)/t_ascend;

    z(i)  = h_target*(10*s^3 - 15*s^4 + 6*s^5);
    vz(i) = h_target*(30*s^2 - 60*s^3 + 30*s^4)/t_ascend;
    az(i) = h_target*(60*s - 180*s^2 + 120*s^3)/t_ascend^2;

end

k = length(t1);

for i = 1:length(t2)

    idx = k+i;

    z(idx)  = h_target;
    vz(idx) = 0.2*sin(5*tgrid(idx));  % small oscillation
    az(idx) = 0;

end

k = length(t1) + length(t2);

for i = 1:length(t3)

    idx = k+i;

    s = (t3(i)-t3(1))/t_descend;

    z(idx)  = h_target*(1-(10*s^3-15*s^4+6*s^5));
    vz(idx) = -h_target*(30*s^2-60*s^3+30*s^4)/t_descend;
    az(idx) = -h_target*(60*s-180*s^2+120*s^3)/t_descend^2;

end

%% THRUST PROFILE (MASS DEPENDENT)

T = m_profile .* (g + az);
T = min(max(T,Tmin),Tmax);

%% ATTITUDE EXCITATION

roll  = deg2rad(2)*sin(0.17*tgrid + pi/4);
pitch = deg2rad(2)*sin(0.20*tgrid);
yaw   = deg2rad(2)*cos(0.15*tgrid);

dt = tgrid(2)-tgrid(1);

P = gradient(roll,dt);
Q = gradient(pitch,dt);
R = gradient(yaw,dt);

%% QUATERNIONS

q = zeros(N,4);

for i = 1:N

phi = roll(i);
theta = pitch(i);
psi = yaw(i);

cy = cos(psi/2);
sy = sin(psi/2);
cp = cos(theta/2);
sp = sin(theta/2);
cr = cos(phi/2);
sr = sin(phi/2);

q(i,:) = [
cr*cp*cy + sr*sp*sy
sr*cp*cy - cr*sp*sy
cr*sp*cy + sr*cp*sy
cr*cp*sy - sr*sp*cy];

end

%% POSITION

X = zeros(N,1);
Y = zeros(N,1);
Z = -z;

Vx = zeros(N,1);
Vy = zeros(N,1);
Vz = -vz;

%% TVC INPUTS

delta_p = deg2rad(1)*sin(0.3*tgrid);
delta_y = deg2rad(1)*cos(0.25*tgrid);

%% RCS INPUT

F_rcs = 20*sin(0.35*tgrid);

%% OUTPUT ARRAYS

Trajectory = [
X Y Z ...
Vx Vy Vz ...
P Q R ...
q];

unom = [
T ...
delta_p ...
delta_y ...
F_rcs];

%% ======================== PLOTS ========================

%% ALTITUDE / VELOCITY / THRUST

figure('Name','Vertical Motion')

subplot(3,1,1)
plot(tgrid,-Z,'LineWidth',2)
ylabel('Altitude (m)')
title('Altitude')
grid on
xlim([0 tf])

subplot(3,1,2)
plot(tgrid,Vz,'LineWidth',2)
ylabel('Vertical Velocity (m/s)')
grid on
xlim([0 tf])

subplot(3,1,3)
plot(tgrid,T,'LineWidth',2)
ylabel('Thrust (N)')
xlabel('Time (s)')
grid on
xlim([0 tf])

%% CONTROL INPUTS

figure('Name','Control Inputs')

subplot(3,1,1)
plot(tgrid,rad2deg(delta_p),'LineWidth',2)
ylabel('Pitch TVC (deg)')
grid on
xlim([0 tf])

subplot(3,1,2)
plot(tgrid,rad2deg(delta_y),'LineWidth',2)
ylabel('Yaw TVC (deg)')
grid on
xlim([0 tf])

subplot(3,1,3)
plot(tgrid,F_rcs,'LineWidth',2)
ylabel('RCS Force (N)')
xlabel('Time (s)')
grid on
xlim([0 tf])

%% ATTITUDE

figure('Name','Attitude')

subplot(3,1,1)
plot(tgrid,rad2deg(roll),'LineWidth',2)
ylabel('Roll (deg)')
grid on
xlim([0 tf])

subplot(3,1,2)
plot(tgrid,rad2deg(pitch),'LineWidth',2)
ylabel('Pitch (deg)')
grid on
xlim([0 tf])

subplot(3,1,3)
plot(tgrid,rad2deg(yaw),'LineWidth',2)
ylabel('Yaw (deg)')
xlabel('Time (s)')
grid on
xlim([0 tf])

%% BODY RATES

figure('Name','Body Rates')

plot(tgrid,[P Q R]*180/pi,'LineWidth',2)
xlabel('Time (s)')
ylabel('Rate (deg/s)')
legend('P','Q','R')
grid on
xlim([0 tf])

%% MASS

figure('Name','Mass Profile')

plot(tgrid,m_profile,'LineWidth',2)
xlabel('Time (s)')
ylabel('Mass (kg)')
grid on
xlim([0 tf])

%% 3D TRAJECTORY

figure('Name','3D Trajectory')

plot3(X,Y,-Z,'LineWidth',2)
xlabel('X (m)')
ylabel('Y (m)')
zlabel('Altitude (m)')
title('Rocket Flight Path')
grid on
view(35,25)

end