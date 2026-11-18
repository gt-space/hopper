function [cases, info] = touchdownCases(results, opts)
%TOUCHDOWNCASES  landingCase per flight Monte Carlo run, from its last state.
%   cases = touchdownCases(R.results)        R from runMonteCarlo
%
%   The flight model stops at first ground contact, so each run's last
%   logged sample is the touchdown state. From it:
%     sink, lateral  NED velocity (vertical down, horizontal speed)
%     heading        drift direction in body axes (azimuth about x from +y)
%     tilt, toward   body x from vertical, and the body azimuth of the low
%                    side (the foot that is lowest)
%     rate           body rates from the last two attitude samples
%     mass           dry mass + remaining propellant
%     thrust         engine thrust at contact, decaying with ThrustTau
%   Runs that errored or diverged are skipped (info.used is false).
%
%   Euler angles are 3-2-1 (yaw, pitch, roll) in degrees, NED, with pitch
%   90 deg upright; the DCM is used, so the singularity there is harmless.

arguments
    results struct
    opts.Vehicle = []
    opts.ThrustTau double = []            % s; default: fitted per run from the cutoff tail-off
end
VC = opts.Vehicle;
if isempty(VC), VC = vehicleConfig('hopper'); end
n = numel(results);
cases = repmat(landingCase(Vehicle=VC), 1, 0);
info = struct('used', false(1, n), 'sink', nan(1, n), 'lateral', nan(1, n), 'tilt', nan(1, n), ...
    'thrust', nan(1, n), 'tau', nan(1, n), 'mass', nan(1, n));
for k = 1:n
    r = results(k);
    if ~isfield(r, 'status') || ~r.status.success || ~isfield(r, 'timeseries') || isempty(r.timeseries)
        continue
    end
    ts = r.timeseries;
    N = numel(ts.t);
    if N < 3, continue; end
    Cbn = dcm(ts.yaw(N), ts.pitch(N), ts.roll(N));      % body -> NED
    Cp  = dcm(ts.yaw(N-1), ts.pitch(N-1), ts.roll(N-1));
    dt = ts.t(N) - ts.t(N-1);
    W = (Cp.' * Cbn - eye(3)) / max(dt, eps);            % ~ skew(w), body
    w = [W(3,2) - W(2,3); W(1,3) - W(3,1); W(2,1) - W(1,2)] / 2;

    vN = [ts.vx(N); ts.vy(N); ts.vz_raw(N)];
    vB = Cbn.' * [vN(1:2); 0];                            % horizontal velocity, body
    dB = Cbn.' * [0; 0; 1];                               % down, body
    tilt = acosd(min(1, -dB(1)));
    toward = atan2d(dB(3), dB(2));
    heading = atan2d(vB(3), vB(2));
    lateral = hypot(vN(1), vN(2));
    mass = VC.mass.dry + ts.ox_mass(N) + ts.fu_mass(N);

    T = ts.thrust(N);
    tau = opts.ThrustTau;
    if isempty(tau), tau = tailOff(ts); end

    c = landingCase(Vehicle=VC, Name=sprintf('run %d', k), SinkRate=max(vN(3), 0), ...
        Lateral=lateral, Heading=heading, Tilt=tilt, TiltToward=toward, Rate=w, ...
        Mass=mass, Thrust=max(T, 0), ThrustTau=tau);
    cases(end+1) = c; %#ok<AGROW>
    info.used(k) = true;
    info.sink(k) = vN(3);  info.lateral(k) = lateral;  info.tilt(k) = tilt;
    info.thrust(k) = T;  info.tau(k) = tau;  info.mass(k) = mass;
end
end

function C = dcm(yaw, pitch, roll)
% body -> NED for 3-2-1 Euler angles in degrees
cy = cosd(yaw); sy = sind(yaw); cp = cosd(pitch); sp = sind(pitch); cr = cosd(roll); sr = sind(roll);
C = [cp*cy, sr*sp*cy - cr*sy, cr*sp*cy + sr*sy;
     cp*sy, sr*sp*sy + cr*cy, cr*sp*sy - sr*cy;
     -sp,   sr*cp,            cr*cp];
end

function tau = tailOff(ts)
% Time constant of the thrust decay over the last 0.3 s (cutoff to contact);
% 0.05 s if the thrust did not fall.
t = ts.t;  T = ts.thrust;
m = t >= t(end) - 0.3;
[Tpk, i] = max(T(m));
tt = t(m);
if Tpk <= 0 || T(end) >= 0.95 * Tpk || tt(end) - tt(i) < 1e-3
    tau = 0.05;
else
    tau = (tt(end) - tt(i)) / log(Tpk / max(T(end), 1e-3 * Tpk));
    tau = min(max(tau, 0.005), 1);
end
end
