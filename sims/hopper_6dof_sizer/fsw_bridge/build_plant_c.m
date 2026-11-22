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
%
%   Output, in Folder:
%     src/          every .c/.h the plant needs (packNGo, flat), so it
%                   builds without MATLAB installed
%     sources.txt   the .c files to compile, one per line
%     defines.txt   the macros MATLAB's build defined (NAME or NAME=VALUE),
%                   one per line
%     plant.json    what was built: vehicle, variants, step, git commit,
%                   gusts, target, ports
%
%   The interface (hopper_6dof.h): hopper_6dof_initialize(), then per
%   1 ms step hopper_6dof_output() (outputs at t) and hopper_6dof_update()
%   (advance to t + step); hopper_6dof_U.u_fsw[4] = [thrust N, TVC pitch
%   rad, TVC yaw rad, RCS torque N m]; hopper_6dof_Y.x_fsw[13] = [pos NED
%   m, vel NED m/s, body rates rad/s, quaternion]. The touchdown stop sets
%   rtmGetStopRequested(hopper_6dof_M); the model keeps stepping if called.
%
%   Parameters are inlined, so changing a workspace value means a rebuild.

here = fileparts(mfilename('fullpath'));
p = inputParser;
p.addParameter('Variants', variantConfig('VEHICLE', 'ctv', 'CONTROLLER', 'fsw'), @isstruct);
p.addParameter('Gusts', true, @islogical);
p.addParameter('Target', 'arm64', @(s) any(strcmp(s, {'arm64', 'x86_64'})));
p.addParameter('Folder', fullfile(here, 'build', 'plant'), @(s) ischar(s) || isstring(s));
p.parse(varargin{:});
V = p.Results.Variants;
out = char(p.Results.Folder);
A = variantAxes();
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
[in, W] = buildSimInput(m, [], [], V);
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
    'built', char(datetime('now', 'Format', 'yyyy-MM-dd''T''HH:mm:ss')), ...
    'inputs', struct('u_fsw', {{'thrust_N', 'tvc_pitch_rad', 'tvc_yaw_rad', 'rcs_torque_Nm'}}), ...
    'outputs', struct('x_fsw', {{'pos_n', 'pos_e', 'pos_d', 'vel_n', 'vel_e', 'vel_d', ...
        'p', 'q', 'r', 'q0', 'q1', 'q2', 'q3'}}));
fid = fopen(fullfile(out, 'plant.json'), 'w');
fprintf(fid, '%s\n', jsonencode(info, 'PrettyPrint', true));
fclose(fid);

rmdir(work, 's');
codeDir = fullfile(out, 'src');
fprintf('Plant code: %d .c files in %s (%.1f s total)\n', numel(srcs), codeDir, toc(t0));
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
