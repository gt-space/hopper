%% Trajectory settings
zf     = 52;     % m
Tup    = 12;     % s
Thover = 2;      % s
Tdown  = 13;     % s
Tterm  = 0;      % terminal descent extension
dt     = 0.01;   % 100 Hz
Tmin = IN.propulsion.nominal_thrust * IN.propulsion.throttle_range(1);
Tmax = 2313;

z0_term = 0;
v0_term = 0;
a0_term = 0;
term_initialized = false;

t_end = Tup + Thover + Tdown + Tterm;
t = (0:dt:t_end)';

%% Mass profile
m0 = OUT.Vehicle.WetMass;        % initial mass
mf = m0-21;        % final mass

m = m0 + (mf-m0)*(t/t_end);   % linear burn

g = 9.81;

%% Minimum jerk profile
mjerk = @(s) (10*s.^3 - 15*s.^4 + 6*s.^5);

%% Build z(t)
z = zeros(size(t));
T_cmd = zeros(size(t));

for k = 1:length(t)

    tk = t(k);

    if tk < Tup

        % ASCENT
        s = tk/Tup;
        z(k) = zf*mjerk(s);

        % upward acceleration
        acc = zf*(60*s - 180*s^2 + 120*s^3)/Tup^2;

        T_cmd(k) = m(k)*(g + acc);

    elseif tk < Tup + Thover

        % HOVER
        z(k) = zf;
        T_cmd(k) = m(k)*g;

    elseif tk < Tup + Thover + Tdown

        % DESCENT
        s = (tk-(Tup+Thover))/Tdown;
        s = min(max(s,0),1);

        z(k) = zf*(1-mjerk(s));

        acc = -zf*(60*s - 180*s^2 + 120*s^3)/Tdown^2;

        T_cmd(k) = m(k)*(g + acc);

    else

        if ~term_initialized
            z0_term = z(k-1);
            v0_term = zdot(k-1);
            a0_term = zddot(k-1);
            term_initialized = true;
        end

        t0 = Tup + Thover + Tdown;
        tau = tk - t0;

        T = Tterm;

        % Solve 5th-order polynomial:
        % z(t) = a0 + a1*t + a2*t^2 + a3*t^3 + a4*t^4 + a5*t^5

        zf_term = 0;     % final position (m)
        vf_term = -0.1;   % final velocity (m/s)
        af_term = 0;      % final acceleration (smooth thrust)

        a0 = zf_term;
        a1 = vf_term;
        a2 = af_term/2;

        % Solve for a3, a4, a5
         A = [ T^3   T^4    T^5;
          3*T^2 4*T^3  5*T^4;
          6*T   12*T^2 20*T^3 ];

         b = [ z0_term - a0 - a1*T - a2*T^2;
          v0_term - a1 - 2*a2*T;
          a0_term - 2*a2 ];             

        x = A\b;

        a3 = x(1);
        a4 = x(2);
        a5 = x(3);

        % Evaluate trajectory
        z(k) = a0 + a1*tau + a2*tau^2 + a3*tau^3 + a4*tau^4 + a5*tau^5;

        zdot_term  = a1 + 2*a2*tau + 3*a3*tau^2 + 4*a4*tau^3 + 5*a5*tau^4;
        zddot_term = 2*a2 + 6*a3*tau + 12*a4*tau^2 + 20*a5*tau^3;

        T_cmd(k) = m(k)*(g + zddot_term);

    end
    zdot  = gradient(z,dt);
    zddot = gradient(zdot,dt);

end

%% Saturate thrust
%% Compute velocity and acceleration
zdot  = gradient(z,dt);
zddot = gradient(zdot,dt);


%% Apply thrust limits
% Tmin = 890;
% Tmax = 2400;

T_cmd = min(max(T_cmd,Tmin),Tmax);

% LOW PASS FILTER 
tau = 0.2; % seconds (tune 0.1–0.5)

T_filt = zeros(size(T_cmd));
T_filt(1) = T_cmd(1);

for k = 2:length(t)
    alpha = dt / (tau + dt);
    T_filt(k) = T_filt(k-1) + alpha*(T_cmd(k) - T_filt(k-1));
end

% Replace thrust command with filtered version
T_cmd = T_filt;

%% Create Signal Editor dataset
T_ts = timeseries(T_cmd,t);
T_ts.Name = "thrust_ref";

scenario = Simulink.SimulationData.Dataset;
scenario = scenario.addElement(T_ts,"thrust_ref");

z_ts = timeseries(z, t); 
z_ts.Name = "z_ref";
vz_ts = timeseries(zdot,t);
scenario = scenario.addElement(z_ts, "z_ref");



% save("thrust_profile.mat","scenario")



% Plots
figure
plot(t,z,"LineWidth",1.5)
grid on
xlabel("Time (s)")
ylabel("z (m)")
title("Reference Position")

figure
plot(t,zddot,"LineWidth",1.5)
grid on
xlabel("Time (s)")
ylabel("Acceleration (m/s^2)")
title("Vertical Acceleration")

figure
plot(t,zdot,"LineWidth",1.5)
grid on
xlabel("Time (s)")
ylabel("Velocity (m/s)")
title("Vertical Velocity")

figure
plot(t,T_cmd,"LineWidth",1.5)
grid on
xlabel("Time (s)")
ylabel("Thrust (N)")
title("Thrust Profile (limited)")
yline(Tmin,'--r')
yline(Tmax,'--r')