% makes trajectorfor_signal_editor.
% Creates a Signal Editor compatible MAT-file that contains:
%   scenario (Simulink.SimulationData.Dataset) with element "z_ref" (timeseries)



% Trajectory settings
zf     = 50;     % m
Tup    = 11.5;      % s  aggressive ascent
Thover = 2;      % s
Tdown  = 13.5;      % s  aggressive descent
dt     = 0.01;   % 100 Hz

t_end = Tup + Thover + Tdown;
t = (0:dt:t_end)';

% Minimum-jerk profile
mjerk = @(s) (10*s.^3 - 15*s.^4 + 6*s.^5);

% Build z(t)
z = zeros(size(t));
for k = 1:numel(t)
    tk = t(k);

    if tk < Tup
        s = tk / Tup;
        z(k) = zf * mjerk(s);

    elseif tk < Tup + Thover
        z(k) = zf;

    else
        s = (tk - (Tup + Thover)) / Tdown;
        s = min(max(s,0),1);
        z(k) = zf * (1 - mjerk(s));
    end
end

% Create Signal Editor "scenario" dataset
z_ts = timeseries(z, t);
z_ts.Name = "z_ref";

scenario = Simulink.SimulationData.Dataset;
scenario = scenario.addElement(z_ts, "z_ref");   % element name inside dataset

% Save MAT file exactly how Signal Editor expects
save("trajectory2.mat","scenario");

% Plot: time vs position
figure;
plot(t, z, "LineWidth", 1.5);
grid on;
xlabel("Time (s)");
ylabel("z_{ref} (m)");
title("Aggressive z_{ref} trajectory (100 Hz)");

disp("Saved Signal Editor compatible file: trajectory2.mat (contains variable: scenario)");

