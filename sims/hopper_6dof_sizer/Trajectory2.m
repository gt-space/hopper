%% Trajectory settings
zf     = 52;     % m
Tup    = 12;     % s
Thover = 2;      % s
Tdown  = 13;     % s
Tterm  = 2.5;    % s, constant-sink terminal segment
z_td   = 0.4;    % m, altitude where the constant-sink segment starts
v_td   = 0.3;    % m/s, touchdown sink rate
dt     = 0.01;   % 100 Hz
% The reference keeps sinking through ground contact so the engine cutoff
% (feet < 0.1 m while descending) always fires and touchdown happens well
% before the end of the finite-horizon gain table, where the attitude gains
% go to zero. Ending at z = 0, v = 0 let the ~0.1 m altitude lag stall the
% vehicle just above the cutoff height in ~7 % of runs (landing transient).

Tmin = IN.propulsion.nominal_thrust * IN.propulsion.throttle_range(1);
Tmax = IN.propulsion.nominal_thrust * IN.propulsion.throttle_range(2);

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

        % DESCENT: quintic from hover to z_td, arriving at sink rate v_td
        s = (tk-(Tup+Thover))/Tdown;

        z(k) = zf + (z_td - zf)*mjerk(s) - v_td*Tdown*(-4*s.^3 + 7*s.^4 - 3*s.^5);

        acc = ((z_td - zf)*(60*s - 180*s^2 + 120*s^3) ...
               - v_td*Tdown*(-24*s + 84*s^2 - 60*s^3))/Tdown^2;

        T_cmd(k) = m(k)*(g + acc);

    else

        % TERMINAL: constant sink through ground contact
        tau = tk - (Tup + Thover + Tdown);
        z(k) = z_td - v_td*tau;
        T_cmd(k) = m(k)*g;

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
% vz_ts = timeseries(zdot,t);
scenario = scenario.addElement(z_ts, "z_ref");



% save("thrust_profile.mat","scenario")



% Plots
% figure
% plot(t,z,"LineWidth",1.5)
% grid on
% xlabel("Time (s)")
% ylabel("z (m)")
% title("Reference Position")
% 
% figure
% plot(t,zddot,"LineWidth",1.5)
% grid on
% xlabel("Time (s)")
% ylabel("Acceleration (m/s^2)")
% title("Vertical Acceleration")
% 
% figure
% plot(t,zdot,"LineWidth",1.5)
% grid on
% xlabel("Time (s)")
% ylabel("Velocity (m/s)")
% title("Vertical Velocity")
% 
% figure
% plot(t,T_cmd,"LineWidth",1.5)
% grid on
% xlabel("Time (s)")
% ylabel("Thrust (N)")
% title("Thrust Profile (limited)")
% yline(Tmin,'--r')
% yline(Tmax,'--r')