function prop = size_propellant(IN, engine, vehicle)

g0 = IN.const.g0;

%% =======================
% Hover requirement
%% =======================
m_guess = vehicle.mass_dry + IN.margins.mass_growth;
T_hover = m_guess * g0;
throttle_hover = T_hover / engine.thrust_nom;

assert(throttle_hover <= max(IN.propulsion.throttle_range), ...
    'Hover throttle exceeds engine capability');

mdot_hover = throttle_hover * engine.mdot_total;
m_hover = mdot_hover * IN.mission.hover_time;

%% =======================
% Ascent Δv (energy-based)
%% =======================
h = IN.mission.target_altitude;
dv_ascent = sqrt(2 * g0 * h) * IN.mission.gravity_loss_factor;

m_ascent = (engine.thrust_nom / (engine.Isp * g0)) * ...
           (dv_ascent / (engine.Isp * g0));

%% =======================
% Descent reserve (soft landing)
%% =======================
dv_descent = 10; % m/s conservative
m_descent = (engine.thrust_nom / (engine.Isp * g0)) * ...
            (dv_descent / (engine.Isp * g0));

%% =======================
% Total propellant
%% =======================
m_prop_total = m_hover + m_ascent + m_descent;

% Split by mixture ratio
OF = IN.propulsion.OF;
m_ox = m_prop_total * OF / (1 + OF);
m_fu = m_prop_total / (1 + OF);

%% =======================
% Pack output
%% =======================
prop = struct();
prop.m_total = m_prop_total;
prop.m_ox = m_ox;
prop.m_fu = m_fu;
prop.mdot_ox = engine.mdot_ox;
prop.mdot_fu = engine.mdot_fu;

end