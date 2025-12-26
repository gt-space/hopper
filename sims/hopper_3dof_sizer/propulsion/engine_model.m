function [thrust, mdot] = engine_model(engine, throttle)

% Enforce throttle bounds
throttle = max(min(throttle,1),0);

% Thrust
thrust = throttle * engine.thrust_nom;

% Mass flow
mdot = throttle * engine.mdot_total;

end
