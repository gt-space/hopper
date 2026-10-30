function leverArmBody = gps_lever_arm(gpsFrame)
% GPS position lever arm: the body-frame vector from the vehicle COM to the
% point whose position the GNSS receiver reports. Used by both the sim
% (gps_model.m) and the filter (gmekf_correct_hopper.m) so they agree.
%
% Inputs:
%
%   gpsFrame    - gpsFrameBus struct with antenna positions r1, r2 (3x1 each)
%
% Outputs:
%
%   leverArmBody - Lever arm in body coordinates (3x1) (m)
%
% NOTE: Dual-antenna (moving base) receivers report the position of the
% primary antenna, antenna 1, plus the baseline to antenna 2. If the chosen
% receiver reports the midpoint instead, use (gpsFrame.r1 + gpsFrame.r2)/2.

    leverArmBody = gpsFrame.r1(:);

end
