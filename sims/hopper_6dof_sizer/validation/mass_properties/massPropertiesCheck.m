% Mass-property model checks
% Run from sims/hopper_6dof_sizer after sim_setup:
%   sim_setup; addpath validation/mass_properties; massPropertiesCheck
% Prints: the mass budget, the derivative tables against finite differences
% of vehicle_mass_props, table integration against direct evaluation, and
% the nominal sim's integrated CG/MOI against vehicle_mass_props at the
% logged tank masses.

mo0 = IN.propulsion.oxidizer_mass; mf0 = IN.propulsion.fuel_mass;
f = @(mo, mf) vehicle_mass_props(MP.DRY, MP.TG, mo, mf);

% ---------- 1. mass budget ----------
fprintf('\n=== MASS BUDGET ===\n');
fprintf('DRY  %s\n', MP.DRY.source);
fprintf('DRY  mass %.3f kg, cg %s m, I diag %s kg m^2\n', MP.DRY.mass, mat2str(MP.DRY.cg', 4), mat2str(diag(MP.DRY.I)', 4));
fprintf('SIM  OUT.Vehicle.WetMass %.3f, DryMass %.3f, Wet-Dry %.3f, propellant loaded %.3f\n', ...
    OUT.Vehicle.WetMass, OUT.Vehicle.DryMass, OUT.Vehicle.WetMass - OUT.Vehicle.DryMass, mo0 + mf0);
[cgF, IF] = f(mo0, mf0); [cgE, IE] = f(0, 0);
fprintf('FULL  cg %.4f  I diag %s\nEMPTY cg %.4f  I diag %s\n', cgF(1), mat2str(diag(IF)', 4), cgE(1), mat2str(diag(IE)', 4));

% ---------- 2. derivative tables vs finite differences ----------
h = 1e-3; E = [];
for fo = [0 0.05 0.25 0.5 0.75 1]
  for ff = [0 0.05 0.25 0.5 0.75 1]
    mo = fo * mo0; mf = ff * mf0;
    dco = interp2(fu_masses, ox_masses, dcg_dMox_table, mf, mo);
    dcf = interp2(fu_masses, ox_masses, dcg_dMfu_table, mf, mo);
    dIo = zeros(1, 3); dIf = dIo; k = [1 5 9];
    for j = 1:3
      dIo(j) = interp2(fu_masses, ox_masses, dI_dMox_table(:, :, k(j)), mf, mo);
      dIf(j) = interp2(fu_masses, ox_masses, dI_dMfu_table(:, :, k(j)), mf, mo);
    end
    a = max(mo - h, 0); b = mo + h; [c1, I1] = f(a, mf); [c2, I2] = f(b, mf);
    fd_co = (c2(1) - c1(1)) / (b - a); fd_Io = diag(I2 - I1)' / (b - a);
    a = max(mf - h, 0); b = mf + h; [c1, I1] = f(mo, a); [c2, I2] = f(mo, b);
    fd_cf = (c2(1) - c1(1)) / (b - a); fd_If = diag(I2 - I1)' / (b - a);
    rel = @(x, y) max(abs(x - y) ./ max(abs(y), 1e-6));
    E(end+1, :) = [rel(dco, fd_co), rel(dIo, fd_Io), rel(dcf, fd_cf), rel(dIf, fd_If)]; %#ok<SAGROW>
  end
end
fprintf('\n=== DERIVATIVE TABLES vs FINITE DIFFERENCE (36 fill points incl. empty) ===\n');
fprintf('DER max rel error: dcg/dMox %.2g%%, dI/dMox %.2g%%, dcg/dMfu %.2g%%, dI/dMfu %.2g%%\n', 100 * max(E, [], 1));

% ---------- 3. table integration vs direct (proportional drain) ----------
s = linspace(1, 0, 400);
cg_i = cgF(1); I_i = diag(IF)';
for k = 2:numel(s)
  dmo = (s(k) - s(k-1)) * mo0; dmf = (s(k) - s(k-1)) * mf0;
  moM = (s(k) + s(k-1)) / 2 * mo0; mfM = (s(k) + s(k-1)) / 2 * mf0;
  cg_i = cg_i + interp2(fu_masses, ox_masses, dcg_dMox_table, mfM, moM) * dmo ...
              + interp2(fu_masses, ox_masses, dcg_dMfu_table, mfM, moM) * dmf;
  for j = 1:3
    kk = [1 5 9];
    I_i(j) = I_i(j) + interp2(fu_masses, ox_masses, dI_dMox_table(:, :, kk(j)), mfM, moM) * dmo ...
                    + interp2(fu_masses, ox_masses, dI_dMfu_table(:, :, kk(j)), mfM, moM) * dmf;
  end
end
fprintf('\n=== TABLE INTEGRATION vs DIRECT (full -> empty) ===\n');
fprintf('INT empty: cg direct %.5f integrated %.5f | I direct %s integrated %s\n', cgE(1), cg_i, mat2str(diag(IE)', 5), mat2str(I_i, 5));

% ---------- 4. nominal sim: integrated outputs vs direct ----------
mdl = 'hopper_6dof_NED_v2'; load_system(mdl);
ph = get_param([mdl '/CG // MOI Script v3'], 'PortHandles');
nm = {'cg_out', 'moi_out', 'didt_out', 'engcg_out'};
for k = 1:4, set_param(ph.Outport(k), 'DataLogging', 'on', 'DataLoggingNameMode', 'Custom', 'DataLoggingName', nm{k}); end
ph1 = get_param([mdl '/Subsystem1'], 'PortHandles');
set_param(ph1.Outport(1), 'DataLogging', 'on', 'DataLoggingNameMode', 'Custom', 'DataLoggingName', 'm_ox');
set_param(ph1.Outport(2), 'DataLogging', 'on', 'DataLoggingNameMode', 'Custom', 'DataLoggingName', 'm_fu');
phs = get_param([mdl '/Saturation1'], 'PortHandles');
set_param(phs.Outport(1), 'DataLogging', 'on', 'DataLoggingNameMode', 'Custom', 'DataLoggingName', 'm_tot');
set_param(mdl, 'SignalLogging', 'on', 'SignalLoggingName', 'logsout');
so = sim(mdl, 'ReturnWorkspaceOutputs', 'on');
g = @(n) getsig(so.logsout, n);
[t, cgo] = g('cg_out'); [~, mo] = g('m_ox'); [~, mf] = g('m_fu'); [~, mt] = g('m_tot'); [~, ec] = g('engcg_out');
[tI, Io] = g('moi_out');
fprintf('\n=== SIM (nominal) ===\n');
fprintf('SIM t_end %.2f s; ox %.3f -> %.3f kg; fu %.3f -> %.3f kg; vehicle mass %.3f -> %.3f kg\n', t(end), mo(1), mo(end), mf(1), mf(end), mt(1), mt(end));
fprintf('SIM vehicle mass vs dry + tanks: max diff %.4f kg\n', max(abs(mt - (OUT.Vehicle.DryMass + mo + mf))));
fprintf('SIM engine CG output at t=0: %s (IN.engine_cg %.4f, IN.mount_cg %.4f)\n', mat2str(ec(1, :), 4), IN.engine_cg, IN.mount_cg);
for p = round(linspace(1, numel(t), 6))
  [cd, Id] = f(mo(p), mf(p));
  pI = find(tI >= t(p), 1); Ii = reshape(Io(pI, :), 3, 3);
  fprintf('SIM t %6.2f  ox %6.3f fu %6.3f | cg_x sim %.4f direct %.4f (diff %+.4f) | I sim %s direct %s\n', ...
      t(p), mo(p), mf(p), cgo(p, 1), cd(1), cgo(p, 1) - cd(1), mat2str(diag(Ii)', 4), mat2str(diag(Id)', 4));
end
close_system(mdl, 0);

function [t, d] = getsig(L, n)
v = L.getElement(n).Values; t = v.Time; d = squeeze(v.Data);
if size(d, 1) ~= numel(t), if ndims(d) == 3, d = reshape(permute(d, [3 1 2]), numel(t), []); else, d = d.'; end, end
end
