%% Trajectory settings (from mission_inputs, so the mission lives in one file)
zf     = IN.mission.target_altitude;  % m
Tup    = IN.mission.ascent_time;      % s
Thover = IN.mission.hover_time;       % s
Tdown  = IN.mission.descent_time;     % s
Tterm  = 0.1;      % terminal descent extension
dt     = 0.01;   % 100 Hz

Tmin = IN.propulsion.nominal_thrust * IN.propulsion.throttle_range(1);
Tmax = IN.propulsion.nominal_thrust * IN.propulsion.throttle_range(2);

z0_term = 0;
v0_term = 0;
a0_term = 0;
term_initialized = false;

t_end = Tup + Thover + Tdown + Tterm;
t = (0:dt:t_end)';

g = 9.81;
% TODO: reconcile with IN.const.g0 (9.80665) - physics change, not done yet

%% Minimum jerk profile
mjerk = @(s) (10*s.^3 - 15*s.^4 + 6*s.^5);

%% Build z(t)
z = zeros(size(t));
acc_cmd = zeros(size(t));   % commanded vertical acceleration; thrust comes later

for k = 1:length(t)

    tk = t(k);

    if tk < Tup

        % ASCENT
        s = tk/Tup;
        z(k) = zf*mjerk(s);

        % upward acceleration
        acc = zf*(60*s - 180*s^2 + 120*s^3)/Tup^2;

        acc_cmd(k) = acc;

    elseif tk < Tup + Thover

        % HOVER
        z(k) = zf;
        acc_cmd(k) = 0;

    elseif tk < Tup + Thover + Tdown

        % DESCENT
        s = (tk-(Tup+Thover))/Tdown;
        s = min(max(s,0),1);

        z(k) = zf*(1-mjerk(s));

        acc = -zf*(60*s - 180*s^2 + 120*s^3)/Tdown^2;

        acc_cmd(k) = acc;

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

        acc_cmd(k) = zddot_term;

    end
    zdot  = gradient(z,dt);
    zddot = gradient(zdot,dt);

end

%% Mass profile
% The propellant burned over the profile, from the commanded thrust and the
% engine mdot lookup, replacing a hardcoded figure. Thrust depends on the mass,
% which depends on the burn, so this iterates; it settles in a few passes.
m0 = OUT.Vehicle.WetMass;        % initial mass

if exist('thrust_data', 'var') && exist('ox_mdot_data', 'var') && exist('fu_mdot_data', 'var')
    mdot_bp  = thrust_data;
    mdot_sum = ox_mdot_data + fu_mdot_data;
else
    mdot_lut = readmatrix('mdot_lookup.xlsx');
    mdot_bp  = mdot_lut(:,1);
    mdot_sum = mdot_lut(:,2) + mdot_lut(:,3);
end

burned = 0;
for iter = 1:20
    m     = m0 - burned*(t/t_end);          % linear burn
    T_try = min(max(m.*(g + acc_cmd), Tmin), Tmax);
    mdot  = interp1(mdot_bp, mdot_sum, ...
        min(max(T_try, mdot_bp(1)), mdot_bp(end)), 'linear');
    burned_new = trapz(t, mdot);
    if abs(burned_new - burned) < 1e-10
        burned = burned_new;
        break
    end
    burned = burned_new;
end

mf = m0 - burned;                % final mass
m  = m0 + (mf-m0)*(t/t_end);     % linear burn
T_cmd = m.*(g + acc_cmd);

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