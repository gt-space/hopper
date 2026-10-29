function [meas, valid, stamp, g] = sensor_gate(meas, valid, stamp, t, faults, g)
% SENSOR_GATE  Apply the fault table's sensor faults to one frame.
%
%   meas, valid, stamp   Sensors block outputs (sensor_layout)
%   faults               fault table (faultScenario)
%   g                    gate state (gate_state): values held at onset
%
% A fault acts while t_on <= t < t_off. Each mode is built to exercise one
% of the GNC raw checks (FDIR CoDR 20-22 / PDR 80): loss -> NaN + invalid,
% saturation -> pinned value, stuck/stale -> held value (and timestamp),
% bias/drift -> offset, noise -> burst, poor_fix -> GPS quality.
%
%#codegen

F = fault_codes();
L = sensor_layout();
for r = 1:size(faults, 1)
    tg = faults(r, 1); md = faults(r, 2);
    t_on = faults(r, 3); t_off = faults(r, 4);
    mag = faults(r, 5); prm = faults(r, 6);
    if tg < F.T_ACCEL || tg > F.T_MAG || t < t_on || t >= t_off
        if t < t_on, g.held(r) = false; end
        continue
    end
    [ch, sens] = targetChannels(tg, prm, L, F);

    % hold the value (and stamp) seen at onset for stuck / stale
    if ~g.held(r)
        g.value(:, r) = meas;
        g.stamp(:, r) = stamp;
        g.held(r) = true;
    end

    switch md
        case F.M_LOSS
            meas(ch) = NaN;
            valid(sens) = false;
        case F.M_BIAS
            meas(ch) = meas(ch) + mag;
        case F.M_DRIFT
            meas(ch) = meas(ch) + mag * (t - t_on);
        case F.M_SATURATION
            meas(ch) = mag;
        case F.M_STUCK
            meas(ch) = g.value(ch, r);
        case F.M_NOISE
            meas(ch) = meas(ch) + mag * sin(2 * pi * 37 * t);   % deterministic 37 Hz burst
        case F.M_STALE
            meas(ch) = g.value(ch, r);
            stamp(sens) = g.stamp(sens, r);
            if tg == F.T_GPS, meas(L.GPS_AGE) = t - g.stamp(L.S_GPS, r); end
        case F.M_POOR_FIX
            meas(L.GPS_FIX) = mag;
            if prm > 0, meas(L.GPS_SATS) = prm; end
        otherwise
    end
end
end

function [ch, sens] = targetChannels(tg, prm, L, F)
% channel mask (meas) and sensor mask (valid / stamp) of a target
ch = false(L.N_MEAS, 1);
sens = false(L.N_SENSORS, 1);
switch tg
    case F.T_ACCEL
        if prm >= 1 && prm <= 3, ch(L.ACCEL(prm)) = true; else, ch(L.ACCEL) = true; end
        sens(L.S_IMU) = true;
    case F.T_GYRO
        if prm >= 1 && prm <= 3, ch(L.GYRO(prm)) = true; else, ch(L.GYRO) = true; end
        sens(L.S_IMU) = true;
    case F.T_BARO
        ch(L.BARO) = true; sens(L.S_BARO) = true;
    case F.T_LIDAR
        if prm >= 1 && prm <= 4
            ch(L.LIDAR(prm)) = true; sens(L.S_LIDAR(prm)) = true;
        else
            ch(L.LIDAR) = true; sens(L.S_LIDAR) = true;
        end
    case F.T_GPS
        ch(L.GPS_POS) = true; sens(L.S_GPS) = true;
    case F.T_MAG
        ch(L.MAG) = true; sens(L.S_MAG) = true;
    otherwise
end
end
