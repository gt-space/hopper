%{

Max Goldstein
AE 6356 Spacecraft Attitude
Problem Set 2
9/13/2026

%}

close all; clear all; clc;

e = [0;0;-1];
e = e./norm(e);

theta = [-720:720,1];

% Euler Rodrigues
q = [sind(theta/2).*e;cosd(theta/2)];

% Alternate Quaternion
qS = -q;

% Rodrigues parameter (Gibbs)
g = e.*tand(theta/2);

% MRP
p = e.*tand(theta/4);

% MRP shadow set
pS = -p./vecnorm(p,2,1).^2;

figure()
plot(theta,q(1,:),'LineWidth',3);
hold on
plot(theta,qS(1,:),'LineWidth',3);
plot(theta,g(1,:),'LineWidth',3);
plot(theta,p(1,:),'LineWidth',3);
plot(theta,pS(1,:),'LineWidth',3);
grid on;
legend('Euler Parameter','Euler Alternate Solution Parameter','Rogrigues Parameter','MRP','MRP Shadow Set');
xlabel('\theta [deg]');
ylabel('Attitude 3rd Component');
title ('\theta vs 1st Component of Attitude Representations')
ylim([-2 2])
%xlim([-180 180])

figure()
plot(theta,q(2,:),'LineWidth',3);
hold on
plot(theta,qS(2,:),'LineWidth',3);
plot(theta,g(2,:),'LineWidth',3);
plot(theta,p(2,:),'LineWidth',3);
plot(theta,pS(2,:),'LineWidth',3);
grid on;
legend('Euler Parameter','Euler Alternate Solution Parameter','Rogrigues Parameter','MRP','MRP Shadow Set');
xlabel('\theta [deg]');
ylabel('Attitude 3rd Component');
title ('\theta vs 2nd Component of Attitude Representations')
ylim([-2 2])
xlim([-180 180])

figure()
plot(theta,q(3,:),'LineWidth',3);
hold on
plot(theta,qS(3,:),'LineWidth',3);
plot(theta,g(3,:),'LineWidth',3);
plot(theta,p(3,:),'LineWidth',3);
plot(theta,pS(3,:),'LineWidth',3);
grid on;
legend('Euler Parameter','Euler Alternate Solution Parameter','Rogrigues Parameter','MRP','MRP Shadow Set');
xlabel('\theta [deg]');
ylabel('Attitude 3rd Component');
title ('\theta vs 3rd Component of Attitude Representations')
ylim([-2 2])
xlim([-180 180])

figure()
plot(theta,p(1,:),'LineWidth',3);
hold on
plot(theta,p(2,:),'LineWidth',3);
plot(theta,p(3,:),'LineWidth',3);
plot(theta,q(1,:),'k--');
plot(theta,q(2,:),'b--');
plot(theta,q(3,:),'r--');
grid on;
legend('p_1','p_2','p_3','q_1','q_2','q_3');
xlabel('\theta [deg]');
ylabel('p component');
title ('\theta vs p')
ylim([-2 2])
%xlim([-180 180])