classdef FlightPhase < Simulink.IntEnumType
    enumeration
        PAD     (1)
        STARTUP (2)
        ASCENT  (3)
        HOVER   (4)
        DESCENT (5)
        LANDED  (6)
    end
end