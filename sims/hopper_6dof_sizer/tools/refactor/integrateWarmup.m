function integrateWarmup()
%INTEGRATEWARMUP  Make the filter warm-up logic work with the FSW plant.
%   The navigation warm-up (PR 50) added a command gate and a liftoff latch
%   to models/hopper_6dof.slx. Two changes for the flight-software plant:
%   - Actuation/Warmup_Control_Gate reads a Digital Clock at the model's
%     fixed step with CONTROLLER=fsw (variant source Clock_Select). Its
%     continuous Clock made the held, discrete command continuous-time, and
%     the generated C then differed from Simulink. Every other controller
%     keeps the Clock: the discrete one moves their first thrust step.
%   - Vehicle_Dynamics: the pad logic (stationary IMU and locked moments
%     until confirmed liftoff, no ground contact after it) ran in every
%     flight. Without a warm-up the IMU read stationary while the vehicle
%     was already moving, and a vehicle on GROUND=strut fell through the
%     pad when it landed. With CONTROLLER=fsw it now runs only with a
%     warm-up configured (controller_start_time > 0). Every other
%     controller keeps it: PR 50's nominal baseline was recorded with it.
%   Rerunning is harmless.

m = 'hopper_6dof';
load_system(m);
cs = getActiveConfigSet(m);
if isa(cs, 'Simulink.ConfigSetRef'), cs = getRefConfigSet(cs); end
changed = false;

% ---- the command gate's clock ----
gate = [m '/Actuation/Warmup_Control_Gate'];
if getSimulinkBlockHandle([gate '/Clock_Select']) < 0
    k = find(variantAxes().CONTROLLER == "fsw");
    assert(isscalar(k), 'integrateWarmup:axis', 'variantAxes has no single fsw choice');
    ph = get_param([gate '/Clock'], 'PortHandles');
    l = get_param(ph.Outport(1), 'Line');
    dst = get_param(l, 'DstPortHandle');
    delete_line(l);
    add_block('simulink/Sources/Digital Clock', [gate '/Clock_FSW'], ...
        'SampleTime', get_param(cs, 'FixedStep'), 'Position', [20 160 50 190]);
    add_block('simulink/Signal Routing/Variant Source', [gate '/Clock_Select'], 'Position', [80 120 110 180]);
    set_param([gate '/Clock_Select'], 'VariantControls', ...
        {sprintf('VARIANT_CONTROLLER ~= %d', k), sprintf('VARIANT_CONTROLLER == %d', k)});
    add_line(gate, 'Clock/1', 'Clock_Select/1', 'autorouting', 'smart');
    add_line(gate, 'Clock_FSW/1', 'Clock_Select/2', 'autorouting', 'smart');
    sp = get_param([gate '/Clock_Select'], 'PortHandles');
    for d = dst(:)', add_line(gate, sp.Outport(1), d, 'autorouting', 'smart'); end
    changed = true;
end

% ---- with fsw, the pad logic needs a warm-up ----
vd = [m '/Vehicle_Dynamics'];
if getSimulinkBlockHandle([vd '/Pad_Logic_On']) < 0
    k = find(variantAxes().CONTROLLER == "fsw");
    assert(isscalar(k), 'integrateWarmup:axis', 'variantAxes has no single fsw choice');
    add_block('simulink/Sources/Constant', [vd '/Pad_Logic_On'], ...
        'Value', sprintf('controller_start_time > 0 || VARIANT_CONTROLLER ~= %d', k), 'OutDataTypeStr', 'boolean', 'Position', [1180 700 1270 720]);
    add_block('simulink/Logic and Bit Operations/Logical Operator', [vd '/Pad_Logic_Off'], ...
        'Operator', 'NOT', 'Position', [1300 695 1330 725]);
    add_block('simulink/Logic and Bit Operations/Logical Operator', [vd '/Off_Pad'], ...
        'Operator', 'OR', 'Inputs', '2', 'Position', [1380 640 1410 680]);
    add_block('simulink/Logic and Bit Operations/Logical Operator', [vd '/Ground_Cutoff'], ...
        'Operator', 'AND', 'Inputs', '2', 'Position', [1380 720 1410 760]);
    add_line(vd, 'Pad_Logic_On/1', 'Pad_Logic_Off/1', 'autorouting', 'smart');
    add_line(vd, 'Airborne_Latch/1', 'Off_Pad/1', 'autorouting', 'smart');
    add_line(vd, 'Pad_Logic_Off/1', 'Off_Pad/2', 'autorouting', 'smart');
    add_line(vd, 'Airborne_Latch/1', 'Ground_Cutoff/1', 'autorouting', 'smart');
    add_line(vd, 'Pad_Logic_On/1', 'Ground_Cutoff/2', 'autorouting', 'smart');
    to = {'Pad_Moment_Lock', 'Off_Pad'; 'Warmup_Accel_Model', 'Off_Pad'; ...
          'Ground_Force_Cutoff', 'Ground_Cutoff'; 'Ground_Moment_Cutoff', 'Ground_Cutoff'};
    for r = 1:size(to, 1)
        ph = get_param([vd '/' to{r, 1}], 'PortHandles');
        src = get_param(get_param(get_param(ph.Inport(2), 'Line'), 'SrcPortHandle'), 'Parent');
        assert(endsWith(src, '/Airborne_Latch'), 'integrateWarmup:wiring', '%s is driven by %s', to{r, 1}, src);
        delete_line(vd, 'Airborne_Latch/1', [to{r, 1} '/2']);
        add_line(vd, [to{r, 2} '/1'], [to{r, 1} '/2'], 'autorouting', 'smart');
    end
    changed = true;
end

if changed
    save_system(m);
    fprintf('%s: fsw: discrete warm-up gate, pad logic only with a warm-up\n', m);
else
    fprintf('%s: warm-up already integrated\n', m);
end
close_system(m, 0);
end
