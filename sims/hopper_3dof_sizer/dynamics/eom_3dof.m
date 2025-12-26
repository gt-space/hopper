function xdot = eom_3dof(~, x, vehicle, prop, thrust, gimbal, IN)

% States
z     = x(1);
zdot  = x(2);
theta = x(3);
q     = x(4);   % pitch rate

m = vehicle.mass;
Iyy = vehicle.Iyy;
g = IN.const.g0;

%% Forces
Tz = thrust * cos(gimbal);
W  = m * g;

zddot = (Tz - W) / m;

%% Moments (about CG)
arm = vehicle.engine_z - vehicle.cg_z;
Myy = thrust * sin(gimbal) * arm;

qdot = Myy / Iyy;

%% Kinematics
thetadot = q;

xdot = [
    zdot;
    zddot;
    thetadot;
    qdot
];

end
