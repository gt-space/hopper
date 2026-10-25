function projectStartup()
%PROJECTSTARTUP  Runs when the HopperSim project opens.
%   Creates the gitignored folders that scripts write into. Simulink build
%   caches go to work/ through the project's cache-folder settings.

for d = {'outputs', fullfile('outputs', 'figures'), fullfile('outputs', 'videos'), 'work'}
    p = hopperRoot(d{1});
    if ~isfolder(p), mkdir(p); end
end
end
