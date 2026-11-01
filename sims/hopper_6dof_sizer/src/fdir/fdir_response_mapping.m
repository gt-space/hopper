function [Mask_Target, Degrade_Multiplier, Abort_Land_Active, ...
          Nav_Fallback_Active, Cutoff_Active, Entry_Trigger] = ...
          fdir_response_mapping(active_response)
% FDIR_RESPONSE_MAPPING  Response tier (0-4) to subsystem command signals.
%
% Tiers follow the GNC PDR ladder (fdir_modes): R1 degrade, R2 abort to
% land, R3 nav-degraded abort to land, R4 cutoff and safing. Only R4 cuts
% the engine.
%
%#codegen

% Initialize outputs to safe nominal defaults
Mask_Target = false;
Degrade_Multiplier = 1.0;
Abort_Land_Active = false;
Nav_Fallback_Active = false;
Cutoff_Active = false;
Entry_Trigger = false;

switch active_response
    case uint8(0) % R0 - Nominal
        % Defaults remain active

    case uint8(1) % R1 - Mask the source, degrade
        Mask_Target = true;
        Degrade_Multiplier = 0.8;

    case uint8(2) % R2 - Abort to land
        Mask_Target = true;
        Degrade_Multiplier = 0.5;
        Abort_Land_Active = true;

    case uint8(3) % R3 - Nav-degraded abort to land (no cutoff)
        Mask_Target = true;
        Degrade_Multiplier = 0.5;
        Abort_Land_Active = true;
        Nav_Fallback_Active = true;

    case uint8(4) % R4 - Cutoff and safing
        Mask_Target = true;
        Degrade_Multiplier = 0.0;
        Abort_Land_Active = true;
        Nav_Fallback_Active = true;
        Cutoff_Active = true;
        Entry_Trigger = true;

    otherwise
        % Default safety fallback
        Mask_Target = true;
        Abort_Land_Active = true;
end
end
