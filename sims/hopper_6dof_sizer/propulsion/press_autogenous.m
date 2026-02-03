function press = press_autogenous(press, tank, dt)
% TODO
R = 188.9; % J/kg-K for N2O gas
T = press.temperature;

% Gas volume increases as liquid leaves
V_g = tank.volume_total - tank.volume_liquid;

% Ideal gas
press.pressure = press.mass_gas * R * T / V_g;

end
