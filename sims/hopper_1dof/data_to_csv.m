thrust = table(out.thrust.Time, out.thrust.Data, ...
          'VariableNames', {'Time','Signal'});
writetable(thrust, 'thrust.csv');

ox_mdot = table(out.ox_mdot.Time, out.ox_mdot.Data, ...
          'VariableNames', {'Time','Signal'});
writetable(ox_mdot, 'ox_mdot.csv');

fuel_mdot = table(out.fuel_mdot.Time, out.fuel_mdot.Data, ...
          'VariableNames', {'Time','Signal'});
writetable(fuel_mdot, 'fu_mdot.csv');