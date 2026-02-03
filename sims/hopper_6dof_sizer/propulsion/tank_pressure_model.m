function P_tank = tank_pressure_model(IN, press, Pc)
% TODO
switch IN.press.mode
    case 'copv'
        % Regulated pressure (ideal regulator)
        P_tank = Pc * (1 + IN.propulsion.injector_dp_frac);

    case 'autogenous'
        P_tank = press.pressure;

    otherwise
        error('Unknown pressurization mode');
end

end
