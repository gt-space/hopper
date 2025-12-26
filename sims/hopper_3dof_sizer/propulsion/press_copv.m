function press = press_copv(press, mdot_gas, dt)

R = 296.8; % J/kg-K for GN2
T = press.temperature;

% Update gas mass
press.mass_gas = press.mass_gas - mdot_gas * dt;

% Ideal gas
press.pressure = press.mass_gas * R * T / press.volume;

end
