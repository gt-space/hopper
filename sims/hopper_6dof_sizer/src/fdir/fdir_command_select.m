function u = fdir_command_select(u_nom, u_lo, x, x_ref, mass, mode, K_abort, g)
% FDIR_COMMAND_SELECT  The command the response ladder lets through.
%
%   u_nom    nominal controller command [thrust; gimbal p; gimbal y; rcs]
%   u_lo     the same controller with the low-bandwidth (R1) gains
%   x, x_ref 13-state and the abort reference (fdir_abort_reference)
%   mass     current vehicle mass (kg)
%   mode     .abort .safed .degrade .nav_abort (fdir_mode_manager /
%            fdir_executive)
%   K_abort  4 x 13 fixed-gain regulator (fdirWorkspace: LQR at hover)
%
% R4 cutoff and safing: thrust 0 (throttle valves closed), gimbal centred,
% RCS off. R2/R3: hover thrust plus the fixed-gain regulator on the abort
% reference; R3 (nav-degraded) without horizontal position and velocity
% feedback, which the distrusted estimate supplies, so it descends on
% attitude, rates and height alone. R1: the low-bandwidth gains. Otherwise
% the nominal command
% passes through unchanged, so a run with FDIR on and no fault flies
% exactly as one with FDIR off.
%
%#codegen

if mode.safed
    u = zeros(4, 1);
elseif mode.abort
    q = x(10:13); qr = x_ref(10:13);
    xe = x(:);
    if dot(q, qr) < 0, xe(10:13) = -q; end     % same quaternion hemisphere
    K = K_abort;
    if isfield(mode, 'nav_abort') && mode.nav_abort
        K(:, [1 2 4 5]) = 0;                    % no N/E position or velocity feedback
    end
    u = [mass * g; 0; 0; 0] - K * (xe - x_ref(:));
elseif mode.degrade
    u = u_lo(:);
else
    u = u_nom(:);
end
end
