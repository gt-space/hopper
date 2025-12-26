function prop = prop_state_update(prop, engine, throttle, dt)
% TODO use aidan code
% Get thrust and mdot
[thrust, mdot] = engine_model(engine, throttle);

% Update propellant masses
prop.m_ox = prop.m_ox - engine.mdot_ox * throttle * dt;
prop.m_fu = prop.m_fu - engine.mdot_fu * throttle * dt;

% Prevent negative mass
prop.m_ox = max(prop.m_ox,0);
prop.m_fu = max(prop.m_fu,0);

% Save outputs
prop.thrust = thrust;
prop.mdot_total = mdot;

end
