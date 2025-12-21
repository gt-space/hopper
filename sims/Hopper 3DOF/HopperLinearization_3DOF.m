function   [A_num, B_num] = HopperLinearization_3DOF(x_val, u_val , data_st)
%% STATES & INPUTS
syms u w q theta xe ze 
x = [u; w; q; theta; xe; ze];

syms T phi  
u_in = [T; phi];


% Vehicle Parameters
vehicle.g = 9.81; % Gravity (m/s^2)

vehicle.Tmax = 1556.88; % max thrust [N]
vehicle.Isp = 260; % seconds
vehicle.mr = 3.8; % O/F Mixture Ratio

vehicle.n2o_mass_0 = 40; % Initial Propellant Masses 
vehicle.ipa_mass_0 = 9; 

vehicle.n2o_tank_height = 0.5; % Propellant tank heights
vehicle.ipa_tank_height = 0.85;

vehicle.n2o_tank_bottom = 0.5; % Bottom position of tanks
vehicle.ipa_tank_bottom = 0.5;

vehicle.cg_dry = 2; % Dry CG of vehicle

vehicle.mass_d = 50; % Dry mass (kg)

vehicle.Iyy_dry = 2500; % Dry moment of inertia (kg*m^2)

vehicle.tvc_pos = 0.1; % Vertical distance of TVC to ground
vehicle.initial_actuator_length = 0.17; % Neutral length of actuator
vehicle.actuator_gain = 4.5; % How much the TVC angle changes with the actuator length (rad/m)

vehicle.ure = 0; % Mass flow reference velocities
vehicle.wre = 0;

g       = vehicle.g;
Isp     = vehicle.Isp;
ure     = vehicle.ure;
wre     = vehicle.wre;
tvc_pos = vehicle.tvc_pos;
mdot = -T / (Isp * g);



[Iyy, Iyy_dot, cg, mdot_n2o, mdot_ipa, m] =  hopper_3dof_inertia_calculator(data_st(7), data_st(8), mdot, vehicle);


% --- TVC Calculations ---
lever = cg - tvc_pos;
Tx =  T * sin(phi); % phi = 0 is TVC no gimbaled, + is to right
Tz = -T * cos(phi);
My =  lever * T * sin(phi);

%% Dynamics:
u_dot = Tx/m - q*w - (mdot/m)*(u - ure) - g*sin(theta);
w_dot = Tz/m + q*u - (mdot/m)*(w - wre) + g*cos(theta);
q_dot = (My - Iyy_dot*q)/Iyy;
theta_dot = q;
xe_dot =  u*cos(theta) - w*sin(theta);
ze_dot =  u*sin(theta) + w*cos(theta);

%% FULL STATE DERIVATIVE
f = [u_dot;
     w_dot;
     q_dot;
     theta_dot;
     xe_dot;
     ze_dot];

%% JACOBIANS: A and B
A = jacobian(f, x);      
B = jacobian(f, u_in);  


A_fun = matlabFunction(A,"Vars",{x, u_in});
B_fun = matlabFunction(B,"Vars",{x, u_in});


A_num = A_fun(x_val, u_val);
B_num = B_fun(x_val, u_val);

end
