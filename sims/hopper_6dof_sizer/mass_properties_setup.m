function MP = mass_properties_setup(IN, TANKS, dry_mass)
% Vehicle mass properties for the sim: a dry block plus the two propellant
% columns
%
%   dry_mass - total dry mass [kg]; the main.m override (VEH.mass.dry)
%   IN.cg_factor, IN.inertia_factor (optional, Monte Carlo) scale the dry
%   CG axial position and the dry inertia
%
%   MP.DRY      dry block: mass, cg [3x1], I [3x3 about dry CG], source
%   MP.TG       as-built tank geometry used for the propellant columns
%   MP.cg_init  full-tank axial CG [m]       (Simulink constant cg_init)
%   MP.MoI_init full-tank inertia [3x3]      (Simulink constant MoI_init)
%   MP.ox_masses, MP.fu_masses, MP.dcg_dMox_table, MP.dcg_dMfu_table,
%   MP.dI_dMox_table, MP.dI_dMfu_table       derivative tables for the
%                                            CG // MOI Script v3 lookups

% ---- Tank geometry (as-built) ----
TG.tank_r   = TANKS.singular.radius;
TG.ox_rho   = 1141;   % kg/m^3
TG.fu_rho   = 800;    % kg/m^3
TG.ox_floor = IN.ox_tank_cg - IN.tanks.ox_h_as_built / 2;
TG.fu_floor = IN.fu_tank_cg - IN.tanks.fu_h_as_built / 2;

% ---- Dry block ----
% Interim, until the CAD mass-properties export exists: dry CG and inertia
% of the cg_moi_init component list with empty tanks, with the inertia
% scaled to dry_mass (keeps its radii of gyration). To use CAD, set
% DRY.cg (datum frame) and DRY.I (about the dry CG, body axes) here.
comp = component_args(IN, TANKS);
[cg_empty, I_empty] = cg_moi_init(0, 0, comp{:});
m_comp = component_mass(comp);

% Monte Carlo dispersion of the dry block (1 in nominal runs): dry CG
% axial position about the datum, and dry inertia
cg_factor = 1; inertia_factor = 1;
if isfield(IN, 'cg_factor'),      cg_factor      = IN.cg_factor;      end
if isfield(IN, 'inertia_factor'), inertia_factor = IN.inertia_factor; end

DRY.mass   = dry_mass;
DRY.cg     = [cg_factor * cg_empty; 0; 0];
DRY.I      = inertia_factor * I_empty * (dry_mass / m_comp);
DRY.source = sprintf('cg_moi_init components (%.2f kg) scaled to %.2f kg', m_comp, dry_mass);
if cg_factor ~= 1 || inertia_factor ~= 1
    DRY.source = sprintf('%s; dispersed cg x%.4f, I x%.4f', DRY.source, cg_factor, inertia_factor);
end

if any(abs(DRY.cg(2:3)) > 0)
    warning('mass_properties_setup:offAxisCG', ...
        'Dry CG is off-axis; the sim CG path uses only the axial component.');
end

% ---- Full-tank values ----
m_ox0 = IN.propulsion.oxidizer_mass;
m_fu0 = IN.propulsion.fuel_mass;
[cg0, I0] = vehicle_mass_props(DRY, TG, m_ox0, m_fu0);

MP.DRY      = DRY;
MP.TG       = TG;
MP.cg_init  = cg0(1);
MP.MoI_init = I0;

% ---- Derivative tables (layout matches the Simulink lookups) ----
n = 100;
MP.ox_masses = linspace(0, m_ox0, n);
MP.fu_masses = linspace(0, m_fu0, n);
MP.dcg_dMox_table = zeros(n, n);
MP.dcg_dMfu_table = zeros(n, n);
MP.dI_dMox_table  = zeros(n, n, 9);
MP.dI_dMfu_table  = zeros(n, n, 9);
for i = 1:n
    for j = 1:n
        [~, ~, dco, dcf, dIo, dIf] = vehicle_mass_props(DRY, TG, MP.ox_masses(i), MP.fu_masses(j));
        MP.dcg_dMox_table(i, j)   = dco;
        MP.dcg_dMfu_table(i, j)   = dcf;
        MP.dI_dMox_table(i, j, :) = reshape(dIo, 1, 1, 9);
        MP.dI_dMfu_table(i, j, :) = reshape(dIf, 1, 1, 9);
    end
end
end

function comp = component_args(IN, TANKS)
% cg_moi_init arguments after (ox_mass, fu_mass)
tank_r = TANKS.singular.radius;
comp = {IN.engine_mass, IN.engine_cg, ...
        TANKS.singular.oxidizer.mass, IN.ox_tank_cg, IN.ox_tank_wall_thick, ...
        TANKS.singular.fuel.mass, IN.fu_tank_cg, IN.fu_tank_wall_thick, ...
        IN.structures.structures_mass, IN.structures_cg, ...
        IN.avionics.avi_mass, IN.avi_cg, ...
        IN.fluids_mass, IN.fluids_cg, ...
        IN.structures.payload_mass, IN.payload_cg, ...
        IN.tanks.ox_h_as_built, IN.tanks.fu_h_as_built, tank_r, tank_r, ...
        IN.press.copv.mass + IN.press.copv.gas_mass, IN.press.copv.r, IN.press.copv.t, ...
        IN.press.copv.h, IN.copv_zcg, IN.copv_ycg};
end

function m = component_mass(comp)
% total mass of the cg_moi_init component list with empty tanks
[engine_mass, ox_tank_mass, ox_wall, fu_tank_mass, fu_wall, structures_mass, ...
 avi_mass, fluids_mass, payload_mass, tank_r, copv_mass] = ...
    deal(comp{1}, comp{3}, comp{5}, comp{6}, comp{8}, comp{9}, comp{11}, comp{13}, comp{15}, comp{19}, comp{21});
al_density = 2700;
caps = 2 * al_density * pi * tank_r^2 * (ox_wall + fu_wall);
m = engine_mass + ox_tank_mass + fu_tank_mass + copv_mass + structures_mass + ...
    avi_mass + fluids_mass + payload_mass + caps;
end
