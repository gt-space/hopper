function PROP = size_propellant(IN, engine, vehicle)

g0 = 9.80665;

% Hover requirement
m_guess = vehicle.mass.dry + IN.margins.mass_growth;
T_hover = m_guess * g0;
throttle_hover = T_hover / IN.propulsion.nominal_thrust;

assert(throttle_hover <= max(IN.propulsion.throttle_range), ...
    'Hover throttle exceeds engine capability');

mdot_hover = throttle_hover * engine.mdot_total; % not exact will be able to better approximate with aidans code
m_hover = mdot_hover * IN.mission.hover_time;

% Ascent delta-v
h = IN.mission.target_altitude;
dv_ascent = sqrt(2 * g0 * h) * IN.mission.gravity_loss_factor;

m_ascent = (IN.propulsion.nominal_thrust / (engine.Isp * g0)) * ... % update calc with aidan code
           (dv_ascent / (engine.Isp * g0));

% Descent reserve (soft landing)
dv_descent = 10; % m/s conservative
m_descent = (IN.propulsion.nominal_thrust / (engine.Isp * g0)) * ... % update calc with aidan codd
            (dv_descent / (engine.Isp * g0));

% Total propellant
m_prop_total = m_hover + m_ascent + m_descent;

% Split by mixture ratio
%OF = IN.propulsion.OF;
%m_ox = m_prop_total * OF / (1 + OF);
%m_fu = m_prop_total / (1 + OF);

% Pack output
PROP = struct();
PROP.m_total = m_prop_total;
PROP.m_ox = m_ox;
PROP.m_fu = m_fu;
PROP.mdot_ox = engine.mdot_ox;
PROP.mdot_fu = engine.mdot_fu;

end