function codeDir = build_env_c(varargin)
%BUILD_ENV_C  Generate and compile C code for the FSW-bridge Environment.
%   build_env_c() copies the Environment subsystem out of
%   hopper_6dof_NED_v2_fswBridge into a standalone model (hopper_env),
%   gives it the bridge model's solver settings, and builds it with
%   Embedded Coder (ert.tlc). Everything generated lands in fsw_bridge/build.
%
%   build_env_c('Gusts', false) turns the random gust disturbance off (used
%   by verify_env_c, since generated code draws different random numbers
%   than Simulink). build_env_c('Folder', dir) builds into another folder.
%
%   Returns the folder holding the generated code.
%
%   Needs the usual base workspace (IN, VEH, lookup tables, ...). If it is
%   missing, sim_setup_cached is run first.
%
%   Parameters are inlined for now, so a change to a workspace value needs a
%   rebuild.

here  = fileparts(mfilename('fullpath'));
sizer = fileparts(here);

p = inputParser;
p.addParameter('Gusts', true, @islogical);
p.addParameter('Folder', fullfile(here, 'build'), @(s) ischar(s) || isstring(s));
p.parse(varargin{:});
build = char(p.Results.Folder);
cd(sizer);

if ~evalin('base', 'exist(''IN'', ''var'') && exist(''VEH'', ''var'')')
    evalin('base', 'sim_setup_cached');
end

src = 'hopper_6dof_NED_v2_fswBridge';
tgt = 'hopper_env';
load_system(src);

% Keep generated code and caches out of the sizer folder
Simulink.fileGenControl('set', 'CodeGenFolder', build, 'CacheFolder', build, 'createDir', true);

% Fresh standalone model holding only the Environment contents
t0 = tic;
if bdIsLoaded(tgt)
    error('build_env_c:loaded', '%s is already loaded; restart MATLAB or close it first.', tgt);
end
new_system(tgt);
Simulink.SubSystem.copyContentsToBlockDiagram([src '/Environment'], tgt);

% The root inport loses its size once nothing drives it, so pin it down:
% u_cmd = [thrust; TVC pitch; TVC yaw; RCS]
set_param([tgt '/u_cmd'], 'PortDimensions', '4', 'OutDataTypeStr', 'double');

if ~p.Results.Gusts
    % Gust on/off comes from this generator; never switching on means no gusts
    set_param([tgt '/Bernoulli Binary Generator'], 'ProbabilityOfZero', '1');
end

cs = copy(getActiveConfigSet(src));
cs.Name = 'fswBridgeConfig';
attachConfigSet(tgt, cs, true);
setActiveConfigSet(tgt, 'fswBridgeConfig');

set_param(tgt, 'SystemTargetFile', 'ert.tlc');
set_param(tgt, ...
    'SupportContinuousTime', 'on', ...   % 6DOF block and integrators
    'SupportAbsoluteTime',   'on', ...   % Clock blocks inside
    'SupportNonFinite',      'on', ...
    'GenerateSampleERTMain', 'on', ...   % example main so the build links
    'GenerateReport',        'off', ...
    'LaunchReport',          'off', ...
    'GenCodeOnly',           'off');
save_system(tgt, fullfile(build, [tgt '.slx']));
fprintf('Standalone model ready in %.1f s\n', toc(t0));

t1 = tic;
slbuild(tgt);
fprintf('Code generation + compile took %.1f s\n', toc(t1));

codeDir = fullfile(build, [tgt '_ert_rtw']);
files = [dir(fullfile(codeDir, '*.c')); dir(fullfile(codeDir, '*.h'))];
fprintf('Generated %d .c/.h files in %s\n', numel(files), codeDir);
end
