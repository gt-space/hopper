function writeVehicleConfig(VC, file, varargin)
%WRITEVEHICLECONFIG  Write a vehicle config struct as a config file.
%   writeVehicleConfig(VC, file)   file is a path ending in <name>_config.m;
%                                  the function it defines is <name>_config
%   writeVehicleConfig(VC, file, 'Sized', {'tanks.ox.mass', ...})
%                                  tag those fields "% sizer <date>"
%
%   For drafts (sizeVehicle) and the CAD importer: the file is a flat list
%   of assignments, every number at full precision, so reading it back
%   gives exactly VC. Review it and fold the values you adopt into the
%   hand-kept config (vehicles/<vehicle>/<vehicle>_config.m), which has
%   the units and reasoning in its comments. Pass the config as its file
%   returns it (feval, not vehicleConfig): the derived fields belong to
%   vehicleConfig, not the file.

p = inputParser;
p.addParameter('Sized', {}, @iscell);
p.parse(varargin{:});
[~, fn] = fileparts(file);
if ~endsWith(fn, '_config')
    error('writeVehicleConfig:name', 'The file name must end in _config.m (got %s).', fn);
end

lines = {sprintf('function VC = %s()', fn)
         sprintf('%%%s  Draft vehicle config written by writeVehicleConfig, %s.', upper(fn), char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm')))
         '%   Review before use: see docs/VEHICLE_CONFIG.md for every field and its units.'
         ''};
lines = [lines; assignments('VC', VC, string(p.Results.Sized))];
lines{end+1} = 'end';
fid = fopen(file, 'w');
if fid < 0, error('writeVehicleConfig:open', 'Cannot write %s.', file); end
closer = onCleanup(@() fclose(fid));
fprintf(fid, '%s\n', lines{:});
end

function L = assignments(prefix, s, sized)
L = {};
for f = fieldnames(s)'
    v = s.(f{1});
    name = [prefix '.' f{1}];
    if isstruct(v) && isscalar(v)
        L = [L; assignments(name, v, sized)]; %#ok<AGROW>
        continue
    end
    line = sprintf('%s = %s;', name, literal(v, name));
    path = extractAfter(string(name), "VC.");
    if any(sized == path)
        line = sprintf('%s   %% sizer %s', line, char(datetime('today', 'Format', 'yyyy-MM')));
    end
    L{end+1, 1} = line; %#ok<AGROW>
end
end

function t = literal(v, name)
if ischar(v)
    t = ['''' strrep(v, '''', '''''') ''''];
elseif isstring(v) && isscalar(v)
    t = ['"' char(strrep(v, '"', '""')) '"'];
elseif islogical(v) && isscalar(v)
    t = mat2str(v);
elseif isnumeric(v) && isreal(v) && ismatrix(v)
    t = mat2str(v, 17);
    if ~isa(v, 'double'), t = sprintf('%s(%s)', class(v), t); end
else
    error('writeVehicleConfig:type', 'Cannot write %s (a %s).', name, class(v));
end
end
