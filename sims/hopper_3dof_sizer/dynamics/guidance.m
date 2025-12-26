function cmd = guidance(state, vehicle, engine, IN)

z    = state.z;
zdot = state.z_dot;
g    = IN.const.g0;

%% Phase detection
if z < IN.mission.target_altitude && zdot >= 0
    phase = 'ascent';
elseif abs(z - IN.mission.target_altitude) < 0.5
    phase = 'hover';
else
    phase = 'descent';
end

%% Throttle logic
switch phase
    case 'ascent'
        Treq = min(engine.thrust_nom, ...
                   vehicle.mass * (g + 2)); % modest accel
    case 'hover'
        Treq = vehicle.mass * g;
    case 'descent'
        % PD controller for landing
        kp = 2.0;
        kd = 1.2;
        Treq = vehicle.mass * ...
            (g - kp*(z - IN.mission.landing_alt) - kd*zdot);
end

throttle = Treq / engine.thrust_nom;
throttle = max(min(throttle,1), IN.propulsion.throttle_range(1));

cmd.throttle = throttle;
cmd.gimbal = 0;  % vertical flight for now

end
