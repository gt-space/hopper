function M = fdir_modes()
% FDIR_MODES  Flight modes, phases and response levels (GNC PDR).
%
% Flight modes follow the GNC PDR state machine (slides 10-11) plus the
% three contingency modes the response ladder drives (slide 79). Phases are
% the fault table's phase mask (CoDR slide 15 / PDR slide 78); each mode
% belongs to one phase. Response levels are R0..R4 (PDR slide 79); R0 and
% R1 never change the mode.
%
%#codegen

% ---- Flight modes ----
M.SAFE              = uint8(1);
M.INITIALIZE        = uint8(2);
M.READY             = uint8(3);
M.ARMED             = uint8(4);
M.ENGINE_START      = uint8(5);
M.PAD_CHECKOUT      = uint8(6);
M.ASCENT            = uint8(7);
M.HOVER             = uint8(8);
M.DESCENT           = uint8(9);
M.TERMINAL_DESCENT  = uint8(10);
M.TOUCHDOWN_CONFIRM = uint8(11);
M.LANDED            = uint8(12);
M.ABORT_LAND        = uint8(13);   % R2
M.NAV_ABORT         = uint8(14);   % R3
M.SAFED             = uint8(15);   % R4
M.NUM_MODES         = 15;

% ---- Phases (fault-table phase mask bits: bit k-1 for phase k) ----
M.PAD     = uint8(1);
M.STARTUP = uint8(2);
M.ASCENT_PHASE  = uint8(3);
M.HOVER_PHASE   = uint8(4);
M.DESCENT_PHASE = uint8(5);
M.LANDED_PHASE  = uint8(6);
M.ALL_PHASES    = uint8(63);       % all six bits

% phase of each mode, indexed by mode number
M.MODE_PHASE = uint8([1 1 1 1 2 2 3 4 5 5 6 6 5 5 6]);

% ---- Response levels ----
M.R0 = uint8(0);   % mask
M.R1 = uint8(1);   % degrade
M.R2 = uint8(2);   % abort to land
M.R3 = uint8(3);   % nav-degraded abort
M.R4 = uint8(4);   % cutoff and safing
M.NO_ARM = uint8(5);   % IBIT failure: hold on the pad (not a ladder level)
end
