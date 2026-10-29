function VC = importCadMassProps(VC, file)
%IMPORTCADMASSPROPS  Dry mass, CG and inertia from a CAD mass-properties export.
%   VC = importCadMassProps(VC, file)   fill VC.mass from the export
%
%   Not implemented yet: the CAD mass-properties export doesn't exist, so
%   its format is unknown. When it does:
%     - read total mass, CG and the inertia tensor about the CG
%     - convert the CAD frame to the body frame (body x up the thrust axis
%       from the bottom centre). The CTV's conversion is the model:
%       body = to_body * CAD with to_body = [0 1 0; 0 0 1; 1 0 0] for a
%       CAD frame with Y up (vehicles/ctv/ctv_config.m, VC.cad)
%     - set VC.mass.dry, VC.mass.dry_cg (axial), VC.mass.cg_lateral,
%       VC.mass.I_full and VC.mass.dry_I, and record the export in
%       VC.source
%     - write the result for review with writeVehicleConfig
%   Until then, CAD values are typed into the config by hand.

error('importCadMassProps:notImplemented', ...
    ['No CAD mass-properties export format yet (%s). Enter the CAD values in ' ...
     'vehicles/<vehicle>/<vehicle>_config.m (VC.mass) for now; see the help text.'], file);
end
