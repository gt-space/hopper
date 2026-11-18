function T = landingSweep(opts)
%LANDINGSWEEP  Crush-core configurations x ground surfaces x touchdown cases.
%   T = landingSweep()   runs the sweep (parfor when a pool can start),
%                        saves results/landing_sweep.mat and returns one
%                        table row per run. landingReport(T) summarizes it.
%
%   Configurations: no core; a core inline in each main strut; a crush
%   pad under each foot; each at a range of plateau forces. Cores are
%   given unlimited length, so the stroke a run uses is the stroke the
%   core needs (divide by the 0.75 usable fraction for the core length).
%
%   opts: MainPlateau (N), PadPlateau (N), Surfaces, PeakRatio, Joints,
%         Mass, Parallel, OutFile

arguments
    opts.MainPlateau double = [6 8 10 12 15 18 22 26 30 35 40] * 1e3
    opts.PadPlateau double = [1.5 2 3 4 5 6 8 10 12] * 1e3
    opts.Surfaces string = ["concrete", "soil_firm", "soil_soft"]
    opts.PeakRatio (1,1) double = 1.0
    opts.Joints (1,:) char = 'pinned'
    opts.Mass (1,1) double = 145
    opts.Parallel (1,1) logical = true
    opts.OutFile (1,:) char = fullfile(hopperRoot(), 'results', 'landing_sweep.mat')
end

VC = vehicleConfig('hopper');
cores = {crushCore('none')};
for F = opts.MainPlateau, cores{end+1} = crushCore('main', Plateau=F, PeakRatio=opts.PeakRatio); end %#ok<AGROW>
for F = opts.PadPlateau,  cores{end+1} = crushCore('pad',  Plateau=F, PeakRatio=opts.PeakRatio); end %#ok<AGROW>
cs = landingCases(Mass=opts.Mass, Vehicle=VC);

[ci, si, ki] = ndgrid(1:numel(cores), 1:numel(opts.Surfaces), 1:numel(cs));
nRun = numel(ci);
fprintf('landingSweep: %d configs x %d surfaces x %d cases = %d runs\n', ...
    numel(cores), numel(opts.Surfaces), numel(cs), nRun);

rows = cell(nRun, 1);
surfaces = opts.Surfaces;  joints = opts.Joints;
usePar = opts.Parallel && canPool();
t0 = tic;
if usePar
    parfor r = 1:nRun
        rows{r} = runOne(VC, cores{ci(r)}, surfaces(si(r)), cs(ki(r)), joints);
    end
else
    for r = 1:nRun
        rows{r} = runOne(VC, cores{ci(r)}, surfaces(si(r)), cs(ki(r)), joints);
    end
end
fprintf('landingSweep: done in %.0f s\n', toc(t0));
T = struct2table(vertcat(rows{:}));
meta = struct('date', char(datetime('now')), 'joints', opts.Joints, 'mass', opts.Mass, ...
    'peak_ratio', opts.PeakRatio, 'vehicle_rev', VC.revision);
outDir = fileparts(opts.OutFile);
if ~exist(outDir, 'dir'), mkdir(outDir); end
save(opts.OutFile, 'T', 'meta');
end

function row = runOne(VC, core, surface, cs, joints)
G = landingGear(Vehicle=VC, Core=core, Joints=joints);
R = landingSim(G, landingSurface(char(surface)), cs.case, Log=false);
pk = R.peak;
row.placement = string(core.placement);
row.plateau_kN = core.plateau / 1e3;
row.surface = surface;
row.case = string(cs.case.name);
row.set = string(cs.set);
row.attitude = string(cs.case.attitude);
row.sink = cs.case.sink;
row.lateral = cs.case.lateral;
row.tilt0 = cs.case.tilt;
row.main_comp_kN = max(pk.Nc(1, :)) / 1e3;
row.main_tens_kN = -min(pk.Nt(1, :)) / 1e3;
row.aux_comp_kN = max(pk.Nc(2:3, :), [], 'all') / 1e3;
row.aux_tens_kN = -min(pk.Nt(2:3, :), [], 'all') / 1e3;
row.main_buckling = max(pk.buckling_ratio(1, :));     % peak compression / Euler load
row.aux_buckling = max(pk.buckling_ratio(2:3, :), [], 'all');
row.sigma_max_MPa = max(pk.sig, [], 'all') / 1e6;
row.foot_Fn_lp_kN = max(pk.Fn_lp) / 1e3;
row.foot_Fn_raw_kN = max(pk.Fn) / 1e3;
row.foot_Ft_kN = max(pk.Ft) / 1e3;
row.grf_total_lp_kN = pk.Fn_total_lp / 1e3;
row.acc_lp_g = pk.acc_lp_g;
row.acc_raw_g = pk.acc_g;
row.stroke_core_mm = pk.stroke_core * 1e3;
row.stroke_pad_mm = pk.stroke_pad * 1e3;
row.dent_mm = pk.dent * 1e3;
row.E_absorbed_J = pk.E_core + pk.E_pad;
row.E_in_J = 0.5 * cs.case.mass * (cs.case.sink^2 + cs.case.lateral^2);
row.tilt_max_deg = pk.tilt_deg;
row.tipped = pk.tipped;
end

function ok = canPool()
ok = false;
try
    ok = license('test', 'Distrib_Computing_Toolbox') && ~isempty(ver('parallel'));
    if ok && isempty(gcp('nocreate')), parpool('Processes'); end
catch
    ok = false;
end
end
