function add_sensors()
%ADD_SENSORS  (Re)build Environment/Sensors in the FSW-bridge model.
%   Creates a Sensors subsystem inside hopper_6dof_NED_v2_fswBridge/Environment
%   holding one MATLAB Function block per sensor, with the code taken from
%   fsw_bridge/sensor_blocks/*.m. Adds Environment outports imu, mag, baro,
%   gps and lidar, logs them at the top level (imu_meas, ...), and saves.
%   Running it again replaces the previous Sensors subsystem, so edit the
%   .m files and re-run rather than editing the blocks by hand.
%
%   Parameters come from the SENS struct (sensor_params). The model's
%   InitFcn creates it if it is missing from the base workspace.

here  = fileparts(mfilename('fullpath'));
sizer = fileparts(here);
cd(sizer);
if ~evalin('base', 'exist(''SENS'', ''var'')')
    assignin('base', 'SENS', sensor_params());
end

mdl = 'hopper_6dof_NED_v2_fswBridge';
env = [mdl '/Environment'];
sen = [env '/Sensors'];
load_system(mdl);

% Sensor blocks: name, inputs (by Sensors inport name), sample time
blocks = {
    'imu_model',   {'F_b', 'W_b', 'mass', 'w_b'},             'SENS.imu.Ts'
    'mag_model',   {'C_bn'},                                  'SENS.mag.Ts'
    'baro_model',  {'p_n', 'v_n', 'C_bn', 'thrust'},          'SENS.baro.Ts'
    'gps_model',   {'p_n', 'v_n', 'C_bn', 'w_b', 'thrust'},   'SENS.gps.tick'
    'lidar_model', {'p_n', 'C_bn'},                           'SENS.lidar.Ts'
    };
outNames = {'imu', 'mag', 'baro', 'gps', 'lidar'};

% Sensors inputs and where they come from inside Environment
inputs = {
    'F_b',    'Sum8',              1   % total force, body axes (N)
    'W_b',    'Product',           1   % weight, body axes (N)
    'mass',   'Saturation1',       1   % vehicle mass (kg)
    'w_b',    '6DOF (Quaternion)', 7   % body rates (rad/s)
    'p_n',    '6DOF (Quaternion)', 2   % CG position NED (m)
    'v_n',    '6DOF (Quaternion)', 1   % CG velocity NED (m/s)
    'C_bn',   '6DOF (Quaternion)', 5   % NED -> body DCM
    'thrust', 'Subsystem2',        1   % delivered thrust (N)
    };

%% Remove a previous version
remove_block([mdl '/log_'], outNames, mdl);
remove_block([env '/'], outNames, env);
if getSimulinkBlockHandle(sen) > 0
    delete_attached_lines(sen);
    delete_block(sen);
end

%% Sensors subsystem
p6 = get_param([env '/6DOF (Quaternion)'], 'Position');
add_block('built-in/Subsystem', sen, 'Position', [p6(1), p6(4) + 150, p6(1) + 160, p6(4) + 400]);
set_param(sen, 'Description', ...
    'Measurement models (Notion: Navigation > Measurement Models). Code: fsw_bridge/sensor_blocks, rebuilt by add_sensors.m. Parameters: SENS (sensor_params.m).');

for i = 1:size(inputs, 1)
    add_block('built-in/Inport', [sen '/' inputs{i, 1}], 'Port', num2str(i), ...
        'Position', [30, 40 + 50 * i, 60, 54 + 50 * i]);
end

rt = sfroot;
for b = 1:size(blocks, 1)
    name = blocks{b, 1};
    path = [sen '/' name];
    add_block('simulink/User-Defined Functions/MATLAB Function', path, ...
        'Position', [250, 40 + 90 * b, 400, 100 + 90 * b]);
    chart = rt.find('-isa', 'Stateflow.EMChart', 'Path', path);
    chart.Script = fileread(fullfile(here, 'sensor_blocks', [name '.m']));
    sensData = chart.find('-isa', 'Stateflow.Data', 'Name', 'SENS');
    sensData.Scope = 'Parameter';
    chart.ChartUpdate = 'DISCRETE';   % sample time only applies to discrete update
    chart.SampleTime = blocks{b, 3};

    add_block('built-in/Outport', [sen '/' outNames{b}], 'Port', num2str(b), ...
        'Position', [480, 62 + 90 * b, 510, 76 + 90 * b]);
    add_line(sen, [name '/1'], [outNames{b} '/1'], 'autorouting', 'on');

    % Chart inputs in port order, wired to the Sensors inport of the same name
    ins = chart.find('-isa', 'Stateflow.Data', 'Scope', 'Input');
    [~, order] = sort(arrayfun(@(d) d.Port, ins));
    ins = ins(order);
    assert(isequal({ins.Name}, blocks{b, 2}), '%s inputs are %s', name, strjoin({ins.Name}, ','));
    for k = 1:numel(ins)
        add_line(sen, [ins(k).Name '/1'], sprintf('%s/%d', name, k), 'autorouting', 'on');
    end
end

%% Wire Sensors into Environment
for i = 1:size(inputs, 1)
    add_line(env, sprintf('%s/%d', inputs{i, 2}, inputs{i, 3}), sprintf('Sensors/%d', i), ...
        'autorouting', 'on');
end
ps = get_param(sen, 'Position');
for b = 1:numel(outNames)
    add_block('built-in/Outport', [env '/' outNames{b}], ...
        'Position', [ps(3) + 120, ps(2) + 40 * b, ps(3) + 150, ps(2) + 14 + 40 * b]);
    add_line(env, sprintf('Sensors/%d', b), [outNames{b} '/1'], 'autorouting', 'on');
end

%% Log the measurements at the top level
pe = get_param(env, 'Position');
for b = 1:numel(outNames)
    port = str2double(get_param([env '/' outNames{b}], 'Port'));
    logger = [mdl '/log_' outNames{b}];
    add_block('simulink/Sinks/To Workspace', logger, 'VariableName', [outNames{b} '_meas'], ...
        'SaveFormat', 'Timeseries', ...
        'Position', [pe(3) + 250, pe(4) + 60 * b, pe(3) + 330, pe(4) + 30 + 60 * b]);
    add_line(mdl, sprintf('Environment/%d', port), ['log_' outNames{b} '/1'], 'autorouting', 'on');
end

%% Make sure SENS exists whenever the model runs
init = get_param(mdl, 'InitFcn');
if ~contains(init, 'sensor_params')
    set_param(mdl, 'InitFcn', strtrim(sprintf('%s\n%s', init, ...
        ['if ~exist(''SENS'', ''var''), addpath(fullfile(fileparts(get_param(bdroot, ''FileName'')), ''fsw_bridge'')); ' ...
         'SENS = sensor_params(); end'])));
end

save_system(mdl);
fprintf('Sensors added to %s (%d sensor blocks)\n', env, size(blocks, 1));
end

function remove_block(prefix, names, ~)
for i = 1:numel(names)
    b = [prefix names{i}];
    if getSimulinkBlockHandle(b) > 0
        delete_attached_lines(b);
        delete_block(b);
    end
end
end

function delete_attached_lines(block)
lh = get_param(block, 'LineHandles');
for L = [lh.Inport(:); lh.Outport(:)]'
    if L > 0 && ishandle(L)
        delete_line(L);
    end
end
end
