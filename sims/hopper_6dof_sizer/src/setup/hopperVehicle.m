function v = hopperVehicle(set)
%HOPPERVEHICLE  The vehicle mission_inputs builds: 'hopper' (default) or 'ctv'.
%   v = hopperVehicle()        current vehicle
%   old = hopperVehicle('ctv') select a vehicle, returning the previous one
%
%   The setup scripts (sim_setup -> main -> mission_inputs) start with
%   clearvars, so the choice can't travel through the base workspace.
%   hopperWorkspace sets it for the duration of one setup and puts the
%   old value back; scripts outside it see 'hopper'.

persistent current
if isempty(current), current = 'hopper'; end
v = current;
if nargin > 0
    A = variantAxes();
    if ~any(A.VEHICLE == string(set))
        error('hopperVehicle:unknown', 'Unknown vehicle ''%s''. Vehicles: %s.', set, strjoin(A.VEHICLE, ', '));
    end
    current = char(set);
end
end
