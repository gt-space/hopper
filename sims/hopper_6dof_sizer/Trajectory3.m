clear
clc

%% ================= USER PARAMETERS =================

m = 114;
g = 9.81;

Tmin = 890;
Tmax = 2400;

h_target = 50;

t_ascend  = 10;
t_hover   = 5;
t_descend = 10;

dt = 0.02;

%% ================= TIME GRID =================

t1 = 0:dt:t_ascend;
t2 = t1(end)+dt:dt:t1(end)+t_hover;
t3 = t2(end)+dt:dt:t2(end)+t_descend;

tgrid = [t1 t2 t3]';
N = length(tgrid);

z  = zeros(N,1);
vz = zeros(N,1);
az = zeros(N,1);

%% ================= MINIMUM JERK ASCENT =================

for i = 1:length(t1)

    s = t1(i)/t_ascend;

    z(i)  = h_target*(10*s^3 - 15*s^4 + 6*s^5);
    vz(i) = h_target*(30*s^2 - 60*s^3 + 30*s^4)/t_ascend;
    az(i) = h_target*(60*s - 180*s^2 + 120*s^3)/t_ascend^2;

end

%% ================= HOVER =================

k = length(t1);

for i = 1:length(t2)

    idx = k+i;

    z(idx)  = h_target;
    vz(idx) = 0;
    az(idx) = 0;

end

%% ================= MINIMUM JERK DESCENT =================

k = length(t1) + length(t2);

for i = 1:length(t3)

    idx = k+i;

    s = (t3(i) - t3(1))/t_descend;

    z(idx)  = h_target*(1 - (10*s^3 - 15*s^4 + 6*s^5));
    vz(idx) = -h_target*(30*s^2 - 60*s^3 + 30*s^4)/t_descend;
    az(idx) = -h_target*(60*s - 180*s^2 + 120*s^3)/t_descend^2;

end

%% ================= THRUST =================

T = m*(g + az);

T = min(max(T,Tmin),Tmax);

%% ================= ATTITUDE EXCITATION =================

roll  = deg2rad(2)*sin(0.17*tgrid + pi/4);
pitch = deg2rad(2)*sin(0.20*tgrid);
yaw   = deg2rad(2)*cos(0.15*tgrid);

%% ================= BODY RATES =================

P = gradient(roll,dt);
Q = gradient(pitch,dt);
R = gradient(yaw,dt);

%% ================= QUATERNIONS =================

q = zeros(N,4);

for i = 1:N

    phi   = roll(i);
    theta = pitch(i);
    psi   = yaw(i);

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
        cr*cp*sy - sr*sp*cy
    ];

end

%% ================= POSITION & VELOCITY =================

X = zeros(N,1);
Y = zeros(N,1);
Z = -z;

Vx = zeros(N,1);
Vy = zeros(N,1);
Vz = -vz;

%% ================= TVC INPUTS =================

delta_p = deg2rad(1)*sin(0.3*tgrid);
delta_y = deg2rad(1)*cos(0.25*tgrid);

F_rcs = zeros(N,1);

%% ================= OUTPUT STRUCTURES =================

Trajectory = [
X Y Z ...
Vx Vy Vz ...
P Q R ...
q
];

unom = [
T
delta_p
delta_y
F_rcs
]';

unom = unom';

%% ================= TRAJECTORY PLOTS =================

figure('Name','Trajectory Overview')

subplot(4,1,1)
plot(tgrid,-Z,'LineWidth',1.5)
ylabel('Altitude (m)')
title('Hopper Trajectory')
grid on

subplot(4,1,2)
plot(tgrid,Vz,'LineWidth',1.5)
ylabel('Vertical Velocity')
grid on

subplot(4,1,3)
plot(tgrid,T,'LineWidth',1.5)
hold on
yline(Tmin,'r--')
yline(Tmax,'r--')
ylabel('Thrust (N)')
legend('Thrust','Min','Max')
grid on

subplot(4,1,4)
plot(tgrid,rad2deg(delta_p),'LineWidth',1.5)
hold on
plot(tgrid,rad2deg(delta_y),'LineWidth',1.5)
ylabel('TVC (deg)')
xlabel('Time (s)')
legend('Pitch','Yaw')
grid on

%% ================= ATTITUDE =================

figure('Name','Attitude')

subplot(3,1,1)
plot(tgrid,rad2deg(roll),'LineWidth',1.5)
ylabel('Roll (deg)')
title('Vehicle Attitude')
grid on

subplot(3,1,2)
plot(tgrid,rad2deg(pitch),'LineWidth',1.5)
ylabel('Pitch (deg)')
grid on

subplot(3,1,3)
plot(tgrid,rad2deg(yaw),'LineWidth',1.5)
ylabel('Yaw (deg)')
xlabel('Time (s)')
grid on

%% ================= BODY RATES =================

figure('Name','Body Rates')

plot(tgrid,[P Q R]*180/pi,'LineWidth',1.5)
xlabel('Time (s)')
ylabel('Angular Rate (deg/s)')
legend('P','Q','R')
title('Body Angular Rates')
grid on

%% ================= QUATERNIONS =================

figure('Name','Quaternion States')

plot(tgrid,q,'LineWidth',1.5)
xlabel('Time (s)')
ylabel('Quaternion Value')
legend('q0','q1','q2','q3')
title('Quaternion States')
grid on

%% ================= 3D FLIGHT PATH =================

figure('Name','3D Flight Path')

plot3(X,Y,-Z,'LineWidth',2)
xlabel('X (m)')
ylabel('Y (m)')
zlabel('Altitude (m)')
title('Rocket Flight Path')
grid on
view(35,25)