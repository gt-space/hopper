%% Trajectory settings
zf     = 50;     % m
Tup    = 12;     % s
Thover = 2;      % s
Tdown  = 13;     % s
Tterm  = 4;      % terminal descent extension
dt     = 0.01;   % 100 Hz
Tmin = 890;
Tmax = 2400;

t_end = Tup + Thover + Tdown + Tterm;
t = (0:dt:t_end)';

%% Mass profile
m0 = OUT.Vehicle.WetMass;        % initial mass
mf = m0-14;        % final mass

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

        % TERMINAL DESCENT
        s = (tk-(Tup+Thover+Tdown))/Tterm;
        s = min(max(s,0),1);

        z_terminal = -1;

        z(k) = (1-s)*0 + s*z_terminal;

        % FORCE MINIMUM THRUST
        T_cmd(k) = Tmin;

    end

end

%% Saturate thrust
%% Compute velocity and acceleration
zdot  = gradient(z,dt);
zddot = gradient(zdot,dt);

%% Required thrust
T = m .* (zddot + g);

%% Apply thrust limits
% Tmin = 890;
% Tmax = 2400;

T_cmd = min(max(T_cmd,Tmin),Tmax);

%% Create Signal Editor dataset
T_ts = timeseries(T_cmd,t);
T_ts.Name = "thrust_ref";

scenario = Simulink.SimulationData.Dataset;
scenario = scenario.addElement(T_ts,"thrust_ref");

z_ts = timeseries(z, t); 
z_ts.Name = "z_ref";

scenario = scenario.addElement(z_ts, "z_ref");



% save("thrust_profile.mat","scenario")



%% Plots
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
plot(t,T_cmd,"LineWidth",1.5)
grid on
xlabel("Time (s)")
ylabel("Thrust (N)")
title("Thrust Profile (limited)")
yline(Tmin,'--r')
yline(Tmax,'--r')