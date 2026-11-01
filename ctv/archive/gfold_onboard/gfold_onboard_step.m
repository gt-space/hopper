function [xref, unom, gfold_flag] = gfold_onboard_step(t, x, m_now, ...
                                        xref_stored, unom_stored, tgrid)
% GFOLD_ONBOARD_STEP  Onboard G-FOLD re-planner, called from Simulink.
%
%   [xref, unom, gfold_flag] = gfold_onboard_step(t, x, m_now, ...
%                                   xref_stored, unom_stored, tgrid)
%
% Runs interpreted (it is invoked through coder.extrinsic from the
% "GFOLD Replan" MATLAB Function block, see gfold_replan.m), so it is free to
% hold persistent state and to call a convex solver.
%
% Behaviour
%   - Ascent and the 51 m hover are NOT re-planned: those rows of the stored
%     reference are passed through untouched. The apex is taken as the last
%     row of xref_stored at peak altitude, i.e. the end of the hover.
%   - From the apex onward, at most once every REPLAN_DT seconds of sim time,
%     it solves a fresh fuel-optimal descent from the CURRENT state
%     (r, v, m) to [0 0 landing_alt] at rest via gfold_onboard, and splices
%     that solution over every tgrid row at or after the current time.
%   - If the solve fails, the previous valid tables are returned unchanged
%     and gfold_flag is raised. That is FDIR fault F-23 (onboard guidance
%     solver failure / tracking a stale plan).
%
% OUTPUTS
%   xref        same size as xref_stored (N x 13, NED + quaternion)
%   unom        same size as unom_stored (N x 4, [T delta_p delta_y F_rcs])
%   gfold_flag  0 = tracking a plan solved at the current state
%               1 = solver failed, tracking the last valid plan (F-23)
%               2 = solver failed and no onboard plan has ever succeeded,
%                   tracking the precomputed reference (F-23)

REPLAN_DT = 1.0;        % s of sim time between solver calls

persistent IN_c xref_c unom_c t_last k_apex have_plan flag_c

%% ---- (re)initialise on first call or on a sim restart -------------------
if isempty(t_last) || t < t_last
    IN_c      = fetch_inputs();
    xref_c    = xref_stored;
    unom_c    = unom_stored;
    t_last    = -inf;
    have_plan = false;
    flag_c    = 0;

    % apex = END of the 51 m hover = last row at peak altitude.
    % xref_stored(:,3) is NED z, so altitude = -z.
    alt_ref = -xref_stored(:,3);
    k_apex  = find(alt_ref >= max(alt_ref) - 1e-6, 1, 'last');
end

% default: hand back the tables and fault state we already hold
xref       = xref_c;
unom       = unom_c;
gfold_flag = flag_c;

if isempty(IN_c)
    % no IN available: pass the precomputed reference through, flagged
    flag_c     = 2;
    gfold_flag = 2;
    return
end

%% ---- only re-plan during descent, and only at the re-plan rate ---------
if t < tgrid(k_apex)
    flag_c     = 0;         % still on the precomputed ascent / hover
    gfold_flag = 0;
    return
end

if (t - t_last) < REPLAN_DT
    return                  % hold the current plan between solver calls
end

t_last = t;

%% ---- stop re-planning once the vehicle is essentially down -------------
alt_now = -x(3);
if alt_now <= IN_c.mission.landing_alt + 0.05
    return
end

%% ---- solve from the current state --------------------------------------
PLAN = gfold_onboard(x(1:3), x(4:6), m_now, IN_c);

if ~PLAN.ok
    % FDIR F-23: keep tracking whatever plan we already had.
    if have_plan
        flag_c = 1;     % holding a previously solved onboard plan
    else
        flag_c = 2;     % never solved; holding the precomputed reference
    end
    gfold_flag = flag_c;
    warning('gfold_onboard_step:solverFailed', ...
        't=%.2f s: onboard G-FOLD failed (%s). Holding last plan (F-23).', ...
        t, PLAN.status);
    return
end

%% ---- splice the new plan onto the stored tables ------------------------
[xref_c, unom_c] = splice(xref_c, unom_c, tgrid, t, PLAN, IN_c);

xref       = xref_c;
unom       = unom_c;
have_plan  = true;
flag_c     = 0;
gfold_flag = 0;

end

% ========================================================================
function [xref, unom] = splice(xref, unom, tgrid, t_now, PLAN, IN)
% Overwrite every reference row at or after t_now with the new plan,
% holding the plan's terminal state past its own end time.

rows = find(tgrid >= t_now);
if isempty(rows)
    return
end

tau = min(max(tgrid(rows) - t_now, 0), PLAN.tf);   % clamp -> hold terminal

r  = interp1(PLAN.t.', PLAN.r_ned.', tau, 'linear', 'extrap');   % n x 3
v  = interp1(PLAN.t.', PLAN.v_ned.', tau, 'linear', 'extrap');   % n x 3
Tv = interp1(PLAN.t.', PLAN.T_ned.', tau, 'linear', 'extrap');   % n x 3

% position / velocity / body rates
xref(rows, 1:3) = r;
xref(rows, 4:6) = v;
xref(rows, 7:9) = 0;            % no commanded body rates
% columns 10:13 (attitude quaternion) keep the stored nominal attitude

% thrust vector -> [T delta_p delta_y 0]
Tmax_gim = IN.tvc.max_gimbal_angle;
for i = 1:numel(rows)
    q     = xref(rows(i), 10:13).';
    [Tm, dp, dy] = thrust_to_gimbal(Tv(i,:).', q, Tmax_gim);
    unom(rows(i), 1) = Tm;
    unom(rows(i), 2) = dp;
    unom(rows(i), 3) = dy;
    unom(rows(i), 4) = 0;       % no roll actuator
end

end

% ========================================================================
function [Tm, delta_p, delta_y] = thrust_to_gimbal(T_ned, q, gim_max)
% Invert the thrust model used by HopperLinearization_6DOF_inertial:
%   Fx_T =  T cos(dp) cos(dy)
%   Fy_T = -T cos(dp) sin(dy)
%   Fz_T = -T sin(dp)

Tm = norm(T_ned);
if Tm < eps
    Tm = 0; delta_p = 0; delta_y = 0;
    return
end

nq = norm(q);
if nq < eps
    q = [1;0;0;0];
else
    q = q / nq;
end

F_Tb = Rbody2NED(q).' * T_ned;      % thrust in body axes

delta_p = asin(min(max(-F_Tb(3)/Tm, -1), 1));
delta_y = atan2(-F_Tb(2), F_Tb(1));

delta_p = min(max(delta_p, -gim_max), gim_max);
delta_y = min(max(delta_y, -gim_max), gim_max);

end

function R = Rbody2NED(q)
q0=q(1); q1=q(2); q2=q(3); q3=q(4);
R = [ ...
 q0^2+q1^2-q2^2-q3^2, 2*(q1*q2-q0*q3),   2*(q1*q3+q0*q2);
 2*(q1*q2+q0*q3),   q0^2-q1^2+q2^2-q3^2, 2*(q2*q3-q0*q1);
 2*(q1*q3-q0*q2),   2*(q2*q3+q0*q1),   q0^2-q1^2-q2^2+q3^2 ];
end

% ========================================================================
function IN = fetch_inputs()
% The model's blocks already resolve IN from the base workspace, so read it
% from there rather than re-running mission_inputs and risking a mismatch.
IN = [];
try
    IN = evalin('base', 'IN');
catch
    try
        IN = mission_inputs();
    catch
        warning('gfold_onboard_step:noInputs', ...
            'Could not resolve IN; onboard G-FOLD disabled, using stored reference.');
    end
end
end
