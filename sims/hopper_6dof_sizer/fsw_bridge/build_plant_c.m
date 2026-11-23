function codeDir = build_plant_c(varargin)
%BUILD_PLANT_C  Generate the plant C code straight from the variant harness.
%   build_plant_c() generates Embedded Coder C for models/hopper_6dof.slx
%   with CONTROLLER=fsw and VEHICLE=ctv: the whole harness minus the
%   controller, which the flight software replaces through the root ports
%   u_fsw (in) and x_fsw (out). Nothing is copied or saved: the variant
%   choices and workspace come from buildSimInput, the code-generation
%   settings go on a copied configuration set, and the model is closed
%   without saving.
%
%   Options:
%     'Variants'  variantConfig(...) result. CONTROLLER must be fsw.
%                 Default variantConfig('VEHICLE', 'ctv', 'CONTROLLER', 'fsw').
%     'Gusts'     false turns the random roll disturbance off (verify_plant_c
%                 needs it off: generated code draws different random
%                 numbers than Simulink). Default true.
%     'Target'    'arm64' (default; the CTV flight computer, LP64) or
%                 'x86_64' (desktop Linux / WSL, also LP64).
%     'Folder'    output folder. Default fsw_bridge/build/plant.
%     'SensorParams'  struct merged over the vehicle's sensor parameters
%                 (IN.sensors, read by sim_sensor_params), e.g.
%                 struct('gps', struct('sigma', [1.5; 1.61; 4.43])).
%                 Default none: the hopper defaults and the vehicle config.
%     'Faults'    fault table (faultScenario(...)) baked into the plant; it
%                 acts only with FAULT_INJECTION on in 'Variants' (sensor
%                 faults in the gates before sens_*, actuator faults on
%                 u_fsw). Default faultScenario(): none.
%     'Ground'    'none' or 'strut': sets GROUND in 'Variants'. Default ''
%                 (keep the Variants' choice, none by default). With 'strut'
%                 the vehicle stands and lands on its leg springs
%                 (subsystems/ground/Ground_Contact.slx) and the touchdown
%                 stop is off (tools/refactor/gateTouchdownStop.m): the
%                 plant keeps flying after contact, and the vehicle comes to
%                 rest, so the flight software's own touchdown detection
%                 sees the impact and the rest period. The plant's engine
%                 cutoff (IN.mission.engine_cutoff_alt, 5 mm on the CTV) still
%                 acts and does not latch, and u_fsw thrust below the
%                 throttle floor is raised to it.
%
%   Output, in Folder:
%     src/          every .c/.h the plant needs (packNGo, flat), so it
%                   builds without MATLAB installed
%     sources.txt   the .c files to compile, one per line
%     defines.txt   the macros MATLAB's build defined (NAME or NAME=VALUE),
%                   one per line
%     plant.json    what was built: vehicle, variants, step, git commit,
%                   gusts, target, ports, and the sensor parameters
%                   and layout baked into the sens_* ports, and
%                   cg_height (the pad CG above the base, m), imu_from_cg
%                   (the IMU from that CG, m, body axes), ground and
%                   touchdown_stop
%
%   The interface (hopper_6dof.h): hopper_6dof_initialize(), then per
%   1 ms step hopper_6dof_output() (outputs at t) and hopper_6dof_update()
%   (advance to t + step); hopper_6dof_U.u_fsw[4] = [thrust N, TVC pitch
%   rad, TVC yaw rad, RCS torque N m]; hopper_6dof_Y.x_fsw[13] = [pos NED
%   m, vel NED m/s, body rates rad/s, quaternion]; hopper_6dof_Y.sens_meas[26],
%   sens_valid[9], sens_stamp[9] = the shared Sensors block's output after
%   the fault gates (FAULT_INJECTION on), laid out as src/sensors/sensor_layout.m:
%   valid is 1/0 per sensor, stamp the time (s) of each sensor's latest
%   sample, which changes when a new sample arrives. The touchdown stop
%   (off with 'Ground', 'strut') sets rtmGetStopRequested(hopper_6dof_M);
%   the model keeps stepping if called.
%
%   Parameters are inlined, so changing a workspace value means a rebuild.

here = fileparts(mfilename('fullpath'));
p = inputParser;
p.addParameter('Variants', variantConfig('VEHICLE', 'ctv', 'CONTROLLER', 'fsw'), @isstruct);
p.addParameter('Gusts', true, @islogical);
p.addParameter('Target', 'arm64', @(s) any(strcmp(s, {'arm64', 'x86_64'})));
p.addParameter('Folder', fullfile(here, 'build', 'plant'), @(s) ischar(s) || isstring(s));
p.addParameter('SensorParams', struct(), @isstruct);
p.addParameter('Faults', faultScenario(), @isnumeric);
p.addParameter('Ground', '', @(s) any(strcmp(s, {'', 'none', 'strut'})));
p.parse(varargin{:});
V = p.Results.Variants;
out = char(p.Results.Folder);
A = variantAxes();
if ~isempty(p.Results.Ground)
    V.VARIANT_GROUND = find(A.GROUND == string(p.Results.Ground));
end
assert(A.CONTROLLER(V.VARIANT_CONTROLLER) == "fsw", 'build_plant_c:controller', ...
    'build_plant_c needs CONTROLLER=fsw (the flight software is the controller).');

m = 'hopper_6dof';
if bdIsLoaded(m)
    error('build_plant_c:loaded', '%s is open; close it (without saving) first.', m);
end
work = fullfile(out, 'work');
if isfolder(out), rmdir(out, 's'); end
mkdir(work);
old = Simulink.fileGenControl('getConfig');
restoreGen = onCleanup(@() Simulink.fileGenControl('setConfig', 'config', old));
Simulink.fileGenControl('set', 'CodeGenFolder', work, 'CacheFolder', work);

t0 = tic;
W = hopperWorkspace([], V);
sensOverride = p.Results.SensorParams;
if ~isempty(fieldnames(sensOverride))
    if ~isfield(W.IN, 'sensors'), W.IN.sensors = struct(); end
    W.IN.sensors = mergeStruct(W.IN.sensors, sensOverride);
end
W.faults = p.Results.Faults;
[in, W] = buildSimInput(m, [], W, V);
load_system(m);
closeModel = onCleanup(@() closeDiscarding(m));
in.applyToModel();
if ~p.Results.Gusts
    % gusts switch on only through this generator; never on = no gusts
    gate = find_system(m, 'MatchFilter', @Simulink.match.allVariants, 'LookUnderMasks', 'all', ...
        'FollowLinks', 'on', 'Name', 'Bernoulli Binary Generator');
    for k = 1:numel(gate), set_param(gate{k}, 'ProbabilityOfZero', '1'); end
end

cs = getActiveConfigSet(m);
if isa(cs, 'Simulink.ConfigSetRef'), cs = getRefConfigSet(cs); end
cs = copy(cs);
cs.Name = 'fswPlantCodegen';
attachConfigSet(m, cs, true);
setActiveConfigSet(m, cs.Name);
hw = struct('arm64', 'ARM Compatible->ARM 64-bit (LP64)', 'x86_64', 'Intel->x86-64 (Linux 64)');
set_param(m, 'SystemTargetFile', 'ert.tlc');
set_param(m, ...
    'ProdHWDeviceType',        hw.(p.Results.Target), ...
    'SupportContinuousTime',   'on', ...   % 6DOF block, integrators, transfer functions
    'SupportAbsoluteTime',     'on', ...   % Clock blocks
    'SupportNonFinite',        'on', ...
    'SupportVariableSizeSignals', 'on', ... % the Sensors block's MATLAB code (sensor_frame) needs it; no malloc results
    'CombineOutputUpdateFcns', 'off', ...  % output() = outputs at t, update() = advance
    'GenerateSampleERTMain',   'off', ...  % the flight software owns main
    'MatFileLogging',          'off', ...
    'GenerateReport',          'off', ...
    'LaunchReport',            'off', ...
    'GenCodeOnly',             'on');
fprintf('Harness set up for codegen in %.1f s\n', toc(t0));

t1 = tic;
slbuild(m);
fprintf('Code generated in %.1f s\n', toc(t1));

% Everything the build needs, flat, so the flight computer needs no MATLAB
rtw = fullfile(work, [m '_ert_rtw']);
bi = load(fullfile(rtw, 'buildInfo.mat'));
zipFile = fullfile(out, 'plant_src.zip');
packNGo(bi.buildInfo, 'packType', 'flat', 'fileName', zipFile);
unzip(zipFile, fullfile(out, 'src'));
delete(zipFile);

srcs = dir(fullfile(out, 'src', '*.c'));
srcs = {srcs.name};
srcs = srcs(~ismember(srcs, {'ert_main.c', 'rt_main.c', 'classic_main.c'}));
writelines(string(srcs(:)), fullfile(out, 'sources.txt'));
defs = string(getDefines(bi.buildInfo));
defs = regexprep(defs(:), '^-D', '');          % NAME or NAME=VALUE
writelines(defs, fullfile(out, 'defines.txt'));

info = struct( ...
    'model', m, ...
    'vehicle', A.VEHICLE(V.VARIANT_VEHICLE), ...
    'controller', A.CONTROLLER(V.VARIANT_CONTROLLER), ...   % "fsw": check this, not the index
    'variants', V, ...
    'fixed_step', str2double(get_param(m, 'FixedStep')), ...
    'fsw_period', W.fsw_period, ...
    'gusts', p.Results.Gusts, ...
    'target', p.Results.Target, ...
    'git_commit', gitCommit(here), ...
    'sensors', true, ...
    'fault_injection', A.FAULT_INJECTION(V.VARIANT_FAULT_INJECTION), ...
    'faults', {num2cell(W.faults, 2)}, ...    % faultScenario rows [target mode t_on t_off magnitude param] (fault_codes); inf -> null
    'sensor_overrides', sensOverride, ...
    'sensor_params', W.sensor_params.S, ...   % sim_sensor_params after the vehicle config and overrides
    'sensor_step', W.sensor_params.dt, ...    % the block's step (the IMU period)
    'sensor_layout', sensor_layout(), ...
    'cg_height', W.cg_init, ...               % full-tank CG above the base on the pad (m, body x; nav_frames heightCOM): sensor geometry is from the base, GPS/LiDAR outputs from the CG
    'imu_from_cg', imuFromCg(W), ...          % IMU from that CG (m, sim body: x up the thrust axis); 0 = the accelerometer reads f at the CG
    'ground', A.GROUND(V.VARIANT_GROUND), ...
    'touchdown_stop', A.GROUND(V.VARIANT_GROUND) ~= "strut", ...   % gateTouchdownStop: off on the struts
    'built', char(datetime('now', 'Format', 'yyyy-MM-dd''T''HH:mm:ss')), ...
    'inputs', struct('u_fsw', {{'thrust_N', 'tvc_pitch_rad', 'tvc_yaw_rad', 'rcs_torque_Nm'}}), ...
    'outputs', struct('x_fsw', {{'pos_n', 'pos_e', 'pos_d', 'vel_n', 'vel_e', 'vel_d', ...
        'p', 'q', 'r', 'q0', 'q1', 'q2', 'q3'}}, ...
        'sens_meas', {sensorChannels()}, ...
        'sens_valid', {sensorNames()}, ...
        'sens_stamp', {sensorNames()}));
fid = fopen(fullfile(out, 'plant.json'), 'w');
fprintf(fid, '%s\n', jsonencode(info, 'PrettyPrint', true));
fclose(fid);

rmdir(work, 's');
codeDir = fullfile(out, 'src');
fprintf('Plant code: %d .c files in %s (%.1f s total)\n', numel(srcs), codeDir, toc(t0));
end

function r = imuFromCg(W)
% S.imu.mount (from the base) relative to the pad CG; zeros when unset (at the CG)
mount = W.sensor_params.S.imu.mount(:);
if all(isfinite(mount))
    r = mount - [W.cg_init; 0; 0];
else
    r = zeros(3, 1);
end
end

function n = sensorChannels()
% names of sens_meas's 26 channels, from sensor_layout
L = sensor_layout();
n = cell(1, L.N_MEAS);
groups = {'ACCEL', {'accel_x', 'accel_y', 'accel_z'}; 'GYRO', {'gyro_x', 'gyro_y', 'gyro_z'}; ...
    'BARO', {'baro_alt'}; 'LIDAR', {'lidar_1', 'lidar_2', 'lidar_3', 'lidar_4'}; ...
    'GPS_POS', {'gps_pos_n', 'gps_pos_e', 'gps_pos_d'}; 'GPS_FIX', {'gps_fix'}; ...
    'GPS_SATS', {'gps_sats'}; 'GPS_AGE', {'gps_age'}; 'MAG', {'mag_x', 'mag_y', 'mag_z'}; ...
    'GPS_VEL', {'gps_vel_n', 'gps_vel_e', 'gps_vel_d'}; 'DGPS', {'dgps_n', 'dgps_e', 'dgps_d'}};
for g = 1:size(groups, 1)
    n(L.(groups{g, 1})) = groups{g, 2};
end
assert(all(~cellfun(@isempty, n)), 'build_plant_c:layout', 'sensor_layout has channels this list misses');
end

function n = sensorNames()
% names of sens_valid / sens_stamp's 9 entries, from sensor_layout
L = sensor_layout();
n = cell(1, L.N_SENSORS);
n(L.S_IMU) = {'imu'}; n(L.S_BARO) = {'baro'}; n(L.S_LIDAR) = {'lidar_1', 'lidar_2', 'lidar_3', 'lidar_4'};
n(L.S_GPS) = {'gps'}; n(L.S_MAG) = {'mag'}; n(L.S_DGPS) = {'dgps'};
assert(all(~cellfun(@isempty, n)), 'build_plant_c:layout', 'sensor_layout has sensors this list misses');
end

function a = mergeStruct(a, b)
for f = fieldnames(b)'
    if isstruct(b.(f{1})) && isfield(a, f{1}) && isstruct(a.(f{1}))
        a.(f{1}) = mergeStruct(a.(f{1}), b.(f{1}));
    else
        a.(f{1}) = b.(f{1});
    end
end
end

function closeDiscarding(m)
close_system(m, 0);
Simulink.data.dictionary.closeAll('-discard');
end

function c = gitCommit(here)
[s, r] = system(sprintf('git -C "%s" rev-parse --short HEAD', here));
if s == 0, c = strtrim(r); else, c = 'unknown'; end
[s, r] = system(sprintf('git -C "%s" status --porcelain -- ..', here));
if s == 0 && ~isempty(strtrim(r)), c = [c '-dirty']; end
end
