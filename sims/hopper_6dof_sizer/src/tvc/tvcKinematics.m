function K = tvcKinematics(P, theta_deg, psi_deg)
%TVCKINEMATICS  Actuator lengths over a pitch/yaw grid (HopperTVCActuatorSizer).
%   K = tvcKinematics(P, -8:0.1:8, -8:0.1:8)
%   K.L_P, K.L_Y   pin-to-pin lengths (m), rows theta, columns psi
%   K.dL_P, K.dL_Y stroke from neutral (m)
%   K.dLdth_P ...  partials (m/rad) for rate mapping
%   K.axis         single-axis tables (other axis 0): delta, dL, r (= dL/ddelta)
%
%   Uses the sizer's constrained-roll gimbal rotation (calculatorGimbal),
%   so both actuators see the coupled motion of a 2-axis gimbal.

if nargin < 2, theta_deg = -8:0.1:8; end
if nargin < 3, psi_deg = theta_deg; end
nT = numel(theta_deg);  nP = numel(psi_deg);
K.theta = theta_deg;  K.psi = psi_deg;
K.L_P = zeros(nT, nP);  K.L_Y = zeros(nT, nP);
for i = 1:nT
    for j = 1:nP
        R = gimbalRotation(theta_deg(i), psi_deg(j));
        K.L_P(i, j) = norm(R * P.RE_P - P.TS_P);
        K.L_Y(i, j) = norm(R * P.RE_Y - P.TS_Y);
    end
end
K.dL_P = K.L_P - P.L0;
K.dL_Y = K.L_Y - norm(P.RE_Y - P.TS_Y);
[K.dLdps_P, K.dLdth_P] = gradient(K.L_P, deg2rad(psi_deg), deg2rad(theta_deg));
[K.dLdps_Y, K.dLdth_Y] = gradient(K.L_Y, deg2rad(psi_deg), deg2rad(theta_deg));

% Single-axis (pitch, yaw = 0) table, fine and wide, for the dynamic sim.
d = deg2rad(-20:0.05:20);
dL = zeros(size(d));
for k = 1:numel(d)
    dL(k) = norm(gimbalRotation(rad2deg(d(k)), 0) * P.RE_P - P.TS_P) - P.L0;
end
K.axis.delta = d;  K.axis.dL = dL;  K.axis.r = gradient(dL, d);
end

function R_GIMB = gimbalRotation(theta, psi)
% calculatorGimbal from HopperTVCActuatorSizer.m (roll-constrained gimbal)
if theta == 0 && psi == 0
    phi_C = 0;
else
    sigma = (cosd(psi)*cosd(theta))/sqrt(sind(theta)^2 + ...
        (cosd(psi)*cosd(theta))^2 + (cosd(theta)*sind(psi))^2);
    gamma = sind(theta)^2 + (cosd(theta)*sind(psi))^2;
    phi_C = atan2d(-(cosd(theta)*sind(psi)*sind(theta)*(sigma - 1))/gamma, ...
        sigma - (cosd(theta)^2*sind(psi)^2*(sigma - 1))/gamma);
end
sx = sind(phi_C); sy = sind(theta); sz = sind(psi);
cx = cosd(phi_C); cy = cosd(theta); cz = cosd(psi);
RX = [1 0 0; 0 cx sx; 0 -sx cx];
RY = [cy 0 -sy; 0 1 0; sy 0 cy];
RZ = [cz sz 0; -sz cz 0; 0 0 1];
R_GIMB = (RX*RY*RZ).';
end
