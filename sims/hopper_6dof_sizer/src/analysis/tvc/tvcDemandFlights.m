function D = tvcDemandFlights(nDispersed, outFile)
%TVCDEMANDFLIGHTS  Gimbal angle / rate / acceleration demand from 6DOF flights.
%   D = tvcDemandFlights(10) flies the nominal hopper plus 10 seeded
%   Monte Carlo scenarios (mc_params.json) with the signal log on, and
%   keeps only the TVC channels (act_bus.tvc_pitch / tvc_yaw, the achieved
%   gimbal angles after the actuator model). Saves D to outFile
%   (default results/tvc_demand.mat) for tvcRequirements / tvcTrade.
%
%   D(k).t, D(k).pitch, D(k).yaw   s, rad (achieved gimbal angle)
%   D(k).name                      'nominal' or 'mc_<k>'

if nargin < 1, nDispersed = 10; end
if nargin < 2, outFile = hopperRoot('results', 'tvc_demand.mat'); end

runs = {struct('name', 'nominal', 'scenario', [])};
if nDispersed > 0
    S = generateScenarios('mc_params.json', nDispersed);
    for k = 1:numel(S)
        runs{end+1} = struct('name', sprintf('mc_%d', k), 'scenario', S(k)); %#ok<AGROW>
    end
end

logDir = fullfile(tempdir, 'tvc_demand_logs');
if ~exist(logDir, 'dir'), mkdir(logDir); end
D = struct('name', {}, 't', {}, 'pitch', {}, 'yaw', {}, 'error', {});
for k = 1:numel(runs)
    fprintf('TVC demand flight %d / %d (%s)\n', k, numel(runs), runs{k}.name);
    logFile = fullfile(logDir, [runs{k}.name '.mat']);
    V = variantConfig('LOG', logFile);
    d = struct('name', runs{k}.name, 't', [], 'pitch', [], 'yaw', [], 'error', '');
    try
        if isempty(runs{k}.scenario)
            out = sim(buildSimInput('hopper_6dof', [], [], V)); %#ok<NASGU>
        else
            out = sim(buildSimInput('hopper_6dof', runs{k}.scenario, [], V)); %#ok<NASGU>
        end
        L = load(logFile);
        act = L.logsout.get('act_bus').Values;
        d.t     = act.tvc_pitch.Time;
        d.pitch = squeeze(act.tvc_pitch.Data);
        d.yaw   = interp1(act.tvc_yaw.Time, squeeze(act.tvc_yaw.Data), d.t);
        delete(logFile);
    catch err
        d.error = err.message;
        fprintf(2, '  failed: %s\n', err.message);
    end
    D(end+1) = d; %#ok<AGROW>
end
save(outFile, 'D');
fprintf('saved %s\n', outFile);
end
