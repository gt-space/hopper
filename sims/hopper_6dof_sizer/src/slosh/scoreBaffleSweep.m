function BS = scoreBaffleSweep(BS)
% Apply the limits in BS to the saved sim metrics and pick, per tank, the
% fewest baffles that pass every damping scale and repeat (tiebreak:
% largest pitch margin). If nothing passes, the design with the smallest worst
% normalized exceedance is picked and BS.stage_pass is false.
% Can be re-run on a loaded baffle_sim_sweep.mat after changing limits.

for c = 1:numel(BS.metrics)
    [BS.metrics(c).pass, BS.metrics(c).worst, BS.metrics(c).pitch_worst] = scoreOne(BS, BS.metrics(c));
end
if isfield(BS, 'final') && ~isempty(BS.final)
    [BS.final.pass, BS.final.worst, BS.final.pitch_worst] = scoreOne(BS, BS.final);
end

for s = 1:2
    idx   = find([BS.cand.stage] == s);
    Nb    = arrayfun(@(c) c.design(s).Nb, BS.cand(idx));
    pass  = [BS.metrics(idx).pass];
    worst = [BS.metrics(idx).worst];
    pworst = [BS.metrics(idx).pitch_worst];

    if any(pass)
        key = [~pass(:), Nb(:) .* pass(:), pworst(:)];  % passing, fewest, most pitch margin
    else
        warning('scoreBaffleSweep:noPass', ...
            'No %s design passes; picking the smallest exceedance.', BS.tank(s).name);
        key = [worst(:), Nb(:)];
    end
    [~, ord] = sortrows(key);
    BS.final_design(s) = BS.cand(idx(ord(1))).design(s);
    BS.stage_pass(s)   = any(pass);
end
end

function [pass, worst, pitch_worst] = scoreOne(BS, M)
% worst: largest metric / limit ratio over all checks (<= 1 passes)
% pitch_worst: largest pitch metric / limit ratio
ratios = [max(M.pitch_err_high(:)) / BS.pitch_lim, ...
          max(M.pitch_err_low(:))  / BS.pitch_lim_low, ...
          max(M.alt_err(:))        / BS.alt_lim, ...
          BS.max_alt_min / min(M.max_alt(:))];
if ~all(M.ok(:))
    ratios(end+1) = Inf;
end
worst = max(ratios);
pitch_worst = max(ratios(1:2));
pass  = worst <= 1;
end
