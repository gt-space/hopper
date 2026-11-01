function parameterizeWindForces()
%PARAMETERIZEWINDFORCES  Wind_Forces reads its aero constants from IN.aero.
%   subsystems/environment/Wind_Forces.slx had the hopper's side area,
%   drag coefficient, air density and centre of pressure as literals. They
%   now come from IN.aero (inputs/mission_inputs.m, same values), so another
%   vehicle (inputs/ctv_inputs.m) can set its own. Rerunning is harmless.

f = 'Wind_Forces';
load_system(f);
map = {'Area', 'IN.aero.area'; 'Cd', 'IN.aero.cd'; 'Density (rho)', 'IN.aero.rho'; 'cp', 'IN.aero.cp'};
for k = 1:size(map, 1)
    set_param([f '/' map{k, 1}], 'Value', map{k, 2});
end
save_system(f);
close_system(f);
fprintf('Wind_Forces: %s now read IN.aero\n', strjoin(map(:, 1)', ', '));
end
