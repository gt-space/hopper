function tests = FDIR_Tests()
% FDIR_TESTS  Unit tests for the FDIR sensor-fault subsystem.
%
%   runtests('FDIR_Tests')
%
% Most of the suite drives the fdir_*.m functions frame by frame, which is
% the same code the MATLAB Function blocks call. The last test runs the
% Simulink model end to end so the wiring is covered too.

tests = functiontests(localfunctions);
end

%% ------------------------- fixture -------------------------

function setupOnce(tc)
tc.TestData.here = fileparts(mfilename('fullpath'));
addpath(tc.TestData.here);
[C, fault_config, fault_map, check_config] = fdir_config();
tc.TestData.C = C;
tc.TestData.fault_config = fault_config;
tc.TestData.fault_map = fault_map;
tc.TestData.check_config = check_config;

% Steady nominal readings: IMU at 1 g down, 50 m baro, 2 m lidar, GPS NED.
tc.TestData.nominal = [0 0 -9.81 0 0 0, 50, 2.0, 10 20 -50]';

% Check and sensor offsets into the 16-wide monitor output.
tc.TestData.IMU_LOSS    = 1;    % IMU       check 1
tc.TestData.IMU_DRIFT   = 2;    % IMU       check 2
tc.TestData.BARO_LOSS   = 5;    % Barometer check 1
tc.TestData.BARO_DRIFT  = 6;    % Barometer check 2
tc.TestData.GPS_STUCK   = 16;   % GPS       check 4
tc.TestData.ASCENT = uint8(3);
tc.TestData.LANDED = uint8(5);
end

%% ------------------------- helpers -------------------------

function r = runFrames(tc, meas, ref, opts)
% Drive monitors -> fault_map -> persistence -> arbitration -> response for
% every column of meas. meas and ref are 11 x nFrames.
C = tc.TestData.C;
fc = tc.TestData.fault_config;
fm = tc.TestData.fault_map;
cc = tc.TestData.check_config;
nFrames = size(meas, 2);

if nargin < 4, opts = struct(); end
if ~isfield(opts, 'flight_state'),    opts.flight_state = repmat(tc.TestData.ASCENT, 1, nFrames); end
if ~isfield(opts, 'ground_reset'),    opts.ground_reset = false(1, nFrames); end
if ~isfield(opts, 'ground_override'), opts.ground_override = false(1, nFrames); end
if ~isfield(opts, 'source_mask'),     opts.source_mask = false(4, nFrames); end

sensorNames = {'IMU', 'Barometer', 'Lidar', 'GPS'};
slices = {1:6, 7, 8, 9:11};
mon = cell(1, 4);
for k = 1:4
    mon{k} = fdir_monitor_state(numel(slices{k}));
end
p_state = fdir_persistence_state(C.NUM_CHECKS, C.M_WINDOW_SIZE);
current_response = uint8(0);
feedback = struct('Response_Ack', false, 'Execution_Status', uint8(0));

r.monitor_fails = false(C.NUM_CHECKS, nFrames);
r.monitor_valid = false(4, nFrames);
r.check_state = zeros(C.NUM_CHECKS, nFrames, 'uint8');   % per monitor check
r.fail_counts = zeros(C.NUM_CHECKS, nFrames, 'uint8');   % per monitor check
r.fault_state = zeros(C.NUM_FAULTS, nFrames, 'uint8');   % per catalog row
r.tier = zeros(1, nFrames, 'uint8');
r.cutoff = false(1, nFrames);

for f = 1:nFrames
    checks = false(16, 1);
    for k = 1:4
        [rf, mv, ~, mon{k}] = fdir_sensor_monitor( ...
            meas(slices{k}, f), ref(slices{k}, f), ...
            opts.source_mask(k, f), opts.ground_reset(f), ...
            C.sensors.(sensorNames{k}), mon{k});
        checks((k-1)*4 + (1:4)) = rf;
        r.monitor_valid(k, f) = mv;
    end
    r.monitor_fails(:, f) = checks;

    % Persistence runs over the 16 checks, then folds onto the 102 rows.
    ground_resets = repmat(opts.ground_reset(f), C.NUM_CHECKS, 1);
    [check_states, counts, p_state] = fdir_persistence( ...
        checks, ground_resets, opts.flight_state(f), cc, p_state);

    fault_states = fdir_fold_checks(check_states, fm, C.NUM_FAULTS);

    [tier, current_response] = fdir_arbitration( ...
        fault_states, opts.flight_state(f), opts.ground_override(f), feedback, ...
        fc, current_response);

    [~, ~, ~, ~, cutoff, ~] = fdir_response_mapping(tier);

    r.check_state(:, f) = [check_states.Persistence_State]';
    r.fault_state(:, f) = [fault_states.Persistence_State]';
    r.fail_counts(:, f) = counts;
    r.tier(f) = tier;
    r.cutoff(f) = cutoff;
end
end

function [meas, ref] = nominalFrames(tc, nFrames)
% Nominal readings with a touch of jitter, so the stuck monitor does not
% trip during runs that are not about the stuck check.
base = tc.TestData.nominal;
ref = repmat(base, 1, nFrames);
meas = ref + 1e-9 * repmat((1:nFrames), numel(base), 1);
end

%% ------------------------- monitor checks -------------------------

function testNominalStaysAtR0(tc)
[meas, ref] = nominalFrames(tc, 30);
r = runFrames(tc, meas, ref);

verifyFalse(tc, any(r.monitor_fails(:)), 'Clean data raised a monitor flag.');
verifyEqual(tc, r.tier, zeros(1, 30, 'uint8'), 'Nominal run left R0.');
verifyFalse(tc, any(r.fault_state(:) == 2), 'Nominal run latched a fault.');
verifyTrue(tc, all(r.monitor_valid(:)), 'Unmasked monitors reported invalid.');
end

function testZeroOnPadIsNotLoss(tc)
% z = 0 sitting on the pad is a real reading. The old check flagged any
% zero as loss of data.
[meas, ref] = nominalFrames(tc, 5);
meas(7, :) = 0;   % barometer altitude
meas(8, :) = 0;   % lidar range
ref(7, :) = 0;
ref(8, :) = 0;
r = runFrames(tc, meas, ref);

verifyFalse(tc, any(r.monitor_fails(5, :)), 'Barometer zero flagged as loss.');
verifyFalse(tc, any(r.monitor_fails(9, :)), 'Lidar zero flagged as loss.');
end

function testEdgeOfRangeIsNotSaturated(tc)
% The old saturation check was inverted and fired on every in-range sample.
C = tc.TestData.C;
[meas, ref] = nominalFrames(tc, 3);
meas(7, :) = C.sensors.Barometer.range_hi;
ref(7, :) = C.sensors.Barometer.range_hi;
meas(8, :) = C.sensors.Lidar.range_lo;
ref(8, :) = C.sensors.Lidar.range_lo;
r = runFrames(tc, meas, ref);

verifyFalse(tc, any(r.monitor_fails(7, :)), 'Barometer at range_hi flagged saturated.');
verifyFalse(tc, any(r.monitor_fails(11, :)), 'Lidar at range_lo flagged saturated.');

% One LSB past the limit must flag.
meas(7, :) = C.sensors.Barometer.range_hi + 1;
r2 = runFrames(tc, meas, ref);
verifyTrue(tc, all(r2.monitor_fails(7, :)), 'Out-of-range barometer not flagged.');
end

function testStuckSensorTripsAfterStuckFrames(tc)
% A frozen GPS was the worst case in testing (206 m flyaway), so the stuck
% check must trip on its own counter.
C = tc.TestData.C;
sf = C.sensors.GPS.stuck_frames;
nFrames = sf + 5;
[meas, ref] = nominalFrames(tc, nFrames);
meas(9:11, :) = repmat(tc.TestData.nominal(9:11), 1, nFrames);   % bit-identical
ref(9:11, :) = meas(9:11, :);
r = runFrames(tc, meas, ref);

stuckRow = r.monitor_fails(tc.TestData.GPS_STUCK, :);
verifyFalse(tc, any(stuckRow(1:sf-1)), 'Stuck flagged before stuck_frames.');
verifyTrue(tc, stuckRow(sf), 'Stuck not flagged on frame stuck_frames.');

verifyEqual(tc, r.check_state(tc.TestData.GPS_STUCK, sf), uint8(2), ...
    'Stuck GPS did not latch on the frame it was detected.');

% F-05 GPS/GNSS dropout is the catalog row GPS stuck reports against.
gpsStuckId = double(tc.TestData.fault_map(tc.TestData.GPS_STUCK));
verifyEqual(tc, gpsStuckId, 5, 'GPS stuck is not mapped to F-05.');
verifyEqual(tc, r.fault_state(gpsStuckId, sf), uint8(2), ...
    'Latched GPS stuck did not fold onto its catalog row.');
end

%% ------------------------- catalog -------------------------

function testCatalogDataRegeneratesFromTheCsv(tc)
% fdir_fault_catalog_data.m is generated from fdir_fault_catalog.csv. If
% someone edits the CSV and forgets to run make_fdir_catalog, or hand-edits
% the generated file, this fails.
here = tc.TestData.here;
csvPath = fullfile(here, 'fdir_fault_catalog.csv');
committed = fullfile(here, 'fdir_fault_catalog_data.m');

verifyTrue(tc, isfile(csvPath), 'fdir_fault_catalog.csv is missing.');
verifyTrue(tc, isfile(committed), 'fdir_fault_catalog_data.m is missing.');

tmp = [tempname '.m'];
cleanup = onCleanup(@() delete(tmp));
make_fdir_catalog(csvPath, tmp);

verifyEqual(tc, fileread(tmp), fileread(committed), ...
    'fdir_fault_catalog_data.m is stale - rerun make_fdir_catalog.');
end

function testSeveritiesMatchTheSheet(tc)
% Spot-check GNC (F-01..F-28) and System (F-98..F-102) against
% FaultConditions.xlsx, sheet "Fault List". The array was previously typed
% in text-sort order (F-10, F-100, F-101, F-102, F-11, ...), which shifted
% every row from F-11 on. It now comes from the sheet via the CSV, so all
% 102 rows are in numeric order, but only GNC and System are spot-checked
% in this PR.
fc = tc.TestData.fault_config;
expected = [11 3; 12 1; 15 5; 19 3; 23 5; 26 2; 98 5; 99 5; 102 5];

for k = 1:size(expected, 1)
    id = expected(k, 1);
    verifyEqual(tc, fc(id).Severity, uint8(expected(k, 2)), ...
        sprintf('F-%02d severity does not match the sheet.', id));
end

% F-73 is N/A in the sheet (throttle valve limit switch, not on the
% rocket), so it carries severity 0 and can never escalate.
verifyEqual(tc, fc(73).Severity, uint8(0), 'F-73 is not marked N/A.');
verifyEqual(tc, fc(73).Response_Tier, uint8(0), 'F-73 does not pin to R0.');

% Tier is severity - 1 across the board; the sheet has no separate
% course-of-action column.
for id = 1:numel(fc)
    verifyEqual(tc, fc(id).Response_Tier, uint8(max(0, double(fc(id).Severity) - 1)), ...
        sprintf('F-%02d tier is not severity - 1.', id));
end
end

function testReadableCatalogAgreesWithTheNumbers(tc)
% The human-facing table and the code-generatable arrays come from the same
% CSV, so they must not drift apart.
T = fdir_fault_catalog();
fc = tc.TestData.fault_config;

verifyEqual(tc, height(T), numel(fc), 'Catalog table and fault_config differ in length.');
verifyEqual(tc, T.ID(15), "F-15", 'Catalog table is not in numeric F-ID order.');

for id = [1 2 4 5 6 7 8 15 73 102]
    verifyEqual(tc, uint8(T.Severity(id)), fc(id).Severity, ...
        sprintf('F-%02d severity differs between the table and fault_config.', id));
    verifyEqual(tc, uint8(T.ResponseTier(id)), fc(id).Response_Tier, ...
        sprintf('F-%02d tier differs between the table and fault_config.', id));
end

verifyTrue(tc, strlength(T.CourseOfAction(15)) > 0, 'F-15 has no course of action.');
end

function testSensorChecksMapToSensorRows(tc)
% Only F-01..F-08 are sensor rows. The placeholder used F-01..F-16, which
% reached F-11 (mass estimate), F-13 (controller instability), F-15 (TWR<1)
% and F-16 (position error).
fm = tc.TestData.fault_map;
verifyEqual(tc, double(fm)', [1 2 2 1, 4 4 4 4, 7 8 7 7, 5 6 6 5], ...
    'fault_map does not match the sheet rows.');
verifyTrue(tc, all(double(fm) <= 8), 'A sensor check maps outside F-01..F-08.');

% No sensor check can command cutoff any more: the sensor rows are S3 (R2)
% except F-06, which is S2 (R1).
fc = tc.TestData.fault_config;
tiers = arrayfun(@(id) double(fc(id).Response_Tier), double(fm));
verifyTrue(tc, all(tiers <= 2), 'A sensor check still escalates past R2.');
end

%% ------------------------- persistence -------------------------

function testTransientDriftDoesNotLatch(tc)
% Three failing frames is below the 4-of-10 threshold. The run is long
% enough for the M = 10 window to slide past the transient: the fault can
% only drop back to INACTIVE once the window holds no failures, which is
% frame 15 for a transient on frames 3 to 5.
nFrames = 18;
[meas, ref] = nominalFrames(tc, nFrames);
meas(1, 3:5) = ref(1, 3:5) + 5.0;   % IMU accel X, well past 5 sigma
r = runFrames(tc, meas, ref);

states = r.check_state(tc.TestData.IMU_DRIFT, :);

verifyEqual(tc, states(3:5), uint8([1 1 1]), 'Transient drift did not go PENDING.');
verifyFalse(tc, any(states == 2), 'Transient drift latched.');
verifyEqual(tc, states(15), uint8(0), 'Transient drift did not return to INACTIVE.');
verifyEqual(tc, r.tier, zeros(1, nFrames, 'uint8'), 'Transient drift escalated the tier.');
end

function testPersistentDriftLatchesOnFourthFrame(tc)
[meas, ref] = nominalFrames(tc, 8);
meas(1, :) = ref(1, :) + 5.0;
r = runFrames(tc, meas, ref);

driftCheck = tc.TestData.IMU_DRIFT;
driftId = double(tc.TestData.fault_map(driftCheck));
states = r.check_state(driftCheck, :);

verifyEqual(tc, states(1:3), uint8([1 1 1]), 'Drift latched before the 4th frame.');
verifyEqual(tc, states(4), uint8(2), 'Drift did not latch on the 4th failing frame.');
verifyEqual(tc, r.fail_counts(driftCheck, 4), uint8(4), 'Wrong fail count at latch.');

% F-02 IMU drift or bias beyond estimate error limits -> severity 3 -> R2.
verifyEqual(tc, driftId, 2, 'IMU drift is not mapped to F-02.');
expectedTier = tc.TestData.fault_config(driftId).Response_Tier;
verifyEqual(tc, r.tier(4), expectedTier, 'Latched drift set the wrong response tier.');
verifyEqual(tc, r.tier(1:3), zeros(1, 3, 'uint8'), 'Tier rose before the latch.');
end

function testNaNLatchesOnFirstFrame(tc)
% A NaN reaches the controller on the next frame, so 4-of-10 persistence
% (40 ms) is too slow. N_Threshold = 1 on the loss faults is the fast path.
[meas, ref] = nominalFrames(tc, 4);
meas(1, :) = NaN;
r = runFrames(tc, meas, ref);

lossCheck = tc.TestData.IMU_LOSS;
lossId = double(tc.TestData.fault_map(lossCheck));
verifyEqual(tc, tc.TestData.check_config(lossCheck).N_Threshold, uint8(1), ...
    'Loss check is not on the fast path.');
verifyTrue(tc, r.monitor_fails(lossCheck, 1), 'NaN not seen as loss.');
verifyEqual(tc, r.check_state(lossCheck, 1), uint8(2), ...
    'Loss did not latch on the first frame.');
verifyEqual(tc, r.tier(1), tc.TestData.fault_config(lossId).Response_Tier, ...
    'Loss latched but the tier did not follow.');
end

function testSharedRowKeepsPerCheckThresholds(tc)
% All four barometer checks report against F-04, but they do not share a
% threshold: loss has to latch in one frame while drift still needs 4 of 10.
% This is why persistence runs over the 16 checks, not the 102 rows.
cc = tc.TestData.check_config;
baroLoss = tc.TestData.BARO_LOSS;
baroDrift = tc.TestData.BARO_DRIFT;

verifyEqual(tc, double(tc.TestData.fault_map(baroLoss)), 4, 'Baro loss is not on F-04.');
verifyEqual(tc, double(tc.TestData.fault_map(baroDrift)), 4, 'Baro drift is not on F-04.');
verifyEqual(tc, cc(baroLoss).N_Threshold, uint8(1), 'Baro loss is not N = 1.');
verifyEqual(tc, cc(baroDrift).N_Threshold, uint8(4), 'Baro drift is not N = 4.');

% Barometer reads NaN: loss latches on frame 1.
[meas, ref] = nominalFrames(tc, 6);
meas(7, :) = NaN;
rLoss = runFrames(tc, meas, ref);
verifyEqual(tc, rLoss.check_state(baroLoss, 1), uint8(2), ...
    'Baro loss did not latch on the first frame.');

% Barometer drifts past 5 sigma: drift needs four frames.
[meas2, ref2] = nominalFrames(tc, 8);
meas2(7, :) = ref2(7, :) + 10.0;
rDrift = runFrames(tc, meas2, ref2);
verifyEqual(tc, rDrift.check_state(baroDrift, 1:3), uint8([1 1 1]), ...
    'Baro drift latched before the 4th frame.');
verifyEqual(tc, rDrift.check_state(baroDrift, 4), uint8(2), ...
    'Baro drift did not latch on the 4th failing frame.');

% Either way F-04 ends up latched, since the row takes the worst state.
verifyEqual(tc, rLoss.fault_state(4, 1), uint8(2), 'F-04 not latched by loss.');
verifyEqual(tc, rDrift.fault_state(4, 4), uint8(2), 'F-04 not latched by drift.');
end

function testInfIsAlsoLoss(tc)
[meas, ref] = nominalFrames(tc, 2);
meas(7, :) = Inf;
r = runFrames(tc, meas, ref);
verifyTrue(tc, r.monitor_fails(5, 1), 'Inf barometer not seen as loss.');
end

%% ------------------------- arbitration -------------------------

function testGroundResetWithOverrideReturnsToR0(tc)
nFrames = 12;
[meas, ref] = nominalFrames(tc, nFrames);
meas(1, 1:6) = ref(1, 1:6) + 5.0;   % latch a drift fault first

opts.flight_state = [repmat(tc.TestData.ASCENT, 1, 6), repmat(tc.TestData.LANDED, 1, 6)];
opts.ground_override = [false(1, 6), true(1, 6)];
opts.ground_reset = [false(1, 6), true(1, 6)];
r = runFrames(tc, meas, ref, opts);

verifyTrue(tc, r.tier(6) > 0, 'Fault never escalated, so the reset proves nothing.');
verifyEqual(tc, r.tier(end), uint8(0), 'Ground reset with override did not return to R0.');

verifyEqual(tc, r.check_state(tc.TestData.IMU_DRIFT, end), uint8(0), ...
    'Ground reset did not clear the latched persistence state.');
end

function testGroundResetWithoutOverrideHoldsTier(tc)
nFrames = 12;
[meas, ref] = nominalFrames(tc, nFrames);
meas(1, 1:6) = ref(1, 1:6) + 5.0;

opts.flight_state = [repmat(tc.TestData.ASCENT, 1, 6), repmat(tc.TestData.LANDED, 1, 6)];
opts.ground_override = false(1, nFrames);
r = runFrames(tc, meas, ref, opts);

verifyTrue(tc, r.tier(end) > 0, 'Tier dropped without the ground override.');
end

function testTierNeverDecreasesMidFlight(tc)
% Escalate on a GPS stuck fault (R4), then hand back clean data. The tier
% must hold for the rest of the flight.
C = tc.TestData.C;
sf = C.sensors.GPS.stuck_frames;
nFrames = sf + 20;
[meas, ref] = nominalFrames(tc, nFrames);
meas(9:11, 1:sf) = repmat(tc.TestData.nominal(9:11), 1, sf);
ref(9:11, 1:sf) = meas(9:11, 1:sf);
r = runFrames(tc, meas, ref);

verifyTrue(tc, max(r.tier) > 0, 'Nothing escalated, so the latch proves nothing.');
verifyTrue(tc, all(diff(double(r.tier)) >= 0), 'Response tier decreased mid-flight.');
verifyEqual(tc, r.tier(end), max(r.tier), 'Tier did not hold at its peak.');
end

function testMaskedSourceReportsInvalidAndNoFaults(tc)
[meas, ref] = nominalFrames(tc, 6);
meas(1, :) = NaN;                       % IMU would otherwise flag loss
opts.source_mask = [true(1, 6); false(3, 6)];
r = runFrames(tc, meas, ref, opts);

verifyFalse(tc, any(r.monitor_valid(1, :)), 'Masked IMU reported valid.');
verifyFalse(tc, any(r.monitor_fails(1:4, :), 'all'), 'Masked IMU raised a fault.');
verifyEqual(tc, r.tier, zeros(1, 6, 'uint8'), 'Masked IMU escalated the tier.');
end

%% ------------------------- response mapping -------------------------

function testResponseMappingTiers(tc)
[mask, deg, abortLand, nav, cutoff, entry] = fdir_response_mapping(uint8(0));
verifyFalse(tc, mask);
verifyEqual(tc, deg, 1.0);
verifyFalse(tc, abortLand || nav || cutoff || entry);

% R3 is a nav-degraded abort, not a mission loss: severity 4 is
% "significant damage to one vehicle subsystem". It still has to fly
% itself down, so the engine stays lit.
[mask3, deg3, abort3, nav3, cutoff3, entry3] = fdir_response_mapping(uint8(3));
verifyTrue(tc, mask3, 'R3 did not mask the source.');
verifyEqual(tc, deg3, 0.5, 'R3 did not derate to 0.5.');
verifyTrue(tc, abort3, 'R3 did not command abort to land.');
verifyTrue(tc, nav3, 'R3 did not command nav fallback.');
verifyFalse(tc, cutoff3, 'R3 cut the engine - an abort to land needs thrust.');
verifyFalse(tc, entry3, 'R3 raised the entry trigger.');

[~, deg4, ~, ~, cutoff4, entry4] = fdir_response_mapping(uint8(4));
verifyEqual(tc, deg4, 0.0);
verifyTrue(tc, cutoff4 && entry4, 'R4 did not command cutoff and entry.');
end

function testOnlyR4CutsTheEngine(tc)
% Cutting thrust is reserved for severity 5, complete loss of mission.
for tier = uint8(0):uint8(4)
    [~, ~, ~, ~, cutoff, ~] = fdir_response_mapping(tier);
    if tier == uint8(4)
        verifyTrue(tc, cutoff, 'R4 did not command cutoff.');
    else
        verifyFalse(tc, cutoff, sprintf('R%d commands cutoff; only R4 may.', tier));
    end
end
end

%% ------------------------- end to end -------------------------

function testModelRunsEndToEnd(tc)
% One sim of the real model: nominal for the first half, then a persistent
% IMU drift that has to escalate through the wired-up blocks.
modelName = 'FDIR_SensorFault_Subsystem';
here = tc.TestData.here;
C = tc.TestData.C;

if ~isfile(fullfile(here, [modelName '.slx']))
    build_fdir_model(modelName);
end

evalin('base', 'Init_FDIR;');

nFrames = 20;
dt = C.SAMPLE_RATE_SEC;
t = (0:nFrames-1)' * dt;
[meas, ref] = nominalFrames(tc, nFrames);
meas(1, 11:end) = ref(1, 11:end) + 5.0;   % IMU accel X drift from frame 11

assignin('base', 'Sensor_data',          [t, meas']);
assignin('base', 'fdir_ekf_ref',         [t, ref']);
assignin('base', 'fdir_source_mask',     [t, zeros(nFrames, 4)]);
assignin('base', 'fdir_ground_reset',    [t, zeros(nFrames, 1)]);
assignin('base', 'fdir_flight_state',    [t, double(tc.TestData.ASCENT) * ones(nFrames, 1)]);
assignin('base', 'fdir_ground_override', [t, zeros(nFrames, 1)]);

if ~bdIsLoaded(modelName)
    load_system(modelName);
end
cleanup = onCleanup(@() close_system(modelName, 0));
set_param(modelName, 'StopTime', num2str((nFrames-1) * dt));
simOut = sim(modelName);

abortLand = tsMatrix(simOut.fdir_log_abort_land);
rawFail = tsMatrix(simOut.fdir_log_raw_fail);

verifyEqual(tc, size(rawFail, 2), C.N_SENSORS * C.N_CHECKS, ...
    'Monitor output is not 4 sensors x 4 checks wide.');
verifyFalse(tc, any(abortLand(1:10)), 'Model escalated during the nominal stretch.');
verifyTrue(tc, any(abortLand(11:end)), 'Model did not escalate on a persistent drift.');

% Drift is severity 3 -> R2 -> abort to land, and R2 does not cut the engine.
verifyFalse(tc, any(tsMatrix(simOut.fdir_log_cutoff)), ...
    'A drift fault commanded engine cutoff.');
end

function d = tsMatrix(ts)
% To Workspace logs a column signal as width x 1 x nFrames; give it back as
% nFrames x width.
raw = ts.Data;
if ndims(raw) == 3
    d = reshape(permute(raw, [3 1 2]), size(raw, 3), []);
else
    d = raw;
end
end
