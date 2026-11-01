function projectStartup()
%PROJECTSTARTUP  Runs when the HopperSim project opens.
%   Creates the gitignored folders that scripts write into. Simulink build
%   caches go to work/ through the project's cache-folder settings.
%   If HOPPER_PYTHON is set (CI does this), MATLAB uses that Python, which
%   must have CoolProp installed (prop_system calls it every MC run).

for d = {'outputs', fullfile('outputs', 'figures'), fullfile('outputs', 'videos'), 'work'}
    p = hopperRoot(d{1});
    if ~isfolder(p), mkdir(p); end
end

py = getenv('HOPPER_PYTHON');
if ~isempty(py) && pyenv().Status == "NotLoaded"
    pyenv(Version=py);
end
end
