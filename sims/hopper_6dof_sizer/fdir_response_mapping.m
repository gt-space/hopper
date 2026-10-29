function [Mask_Target, Degrade_Multiplier, Abort_Land_Active, ...
          Nav_Fallback_Active, Cutoff_Active, Entry_Trigger] = ...
          fdir_response_mapping(active_response)
% FDIR_RESPONSE_MAPPING  Response tier (0-4) to subsystem command signals.
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

    case uint8(1) % R1 - Reconfiguration / Minor Degradation
        Mask_Target = true;
        Degrade_Multiplier = 0.8;

    case uint8(2) % R2 - Safe Mode / Abort Landing
        Mask_Target = true;
        Degrade_Multiplier = 0.5;
        Abort_Land_Active = true;

    case uint8(3) % R3 - Nav-Degraded Abort
        % Severity 4 in FaultConditions.xlsx is "significant damage to one
        % vehicle subsystem", not "complete loss of mission" (severity 5).
        % An R3 abort still has to fly itself down, so the engine stays lit
        % and throttle is derated rather than cut. Only R4 cuts.
        Mask_Target = true;
        Degrade_Multiplier = 0.5;
        Abort_Land_Active = true;
        Nav_Fallback_Active = true;

    case uint8(4) % R4 - Emergency Escape / Abort
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
