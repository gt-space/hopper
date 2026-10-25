function [active_response, current_response] = fdir_arbitration( ...
        fault_states, flight_state, ground_override, response_feedback, ...
        fault_config, current_response)
% FDIR_ARBITRATION  Highest requested response tier, latched.
%
% response_feedback carries Response_Feedback_Bus. Nothing in this subsystem
% produces an acknowledgement yet, so the input is accepted and not read;
% it stays in the signature so the port exists once a responder is wired in.
%
%   active_response   current tier, 0 (R0) .. 4 (R4)
%   current_response  latched tier, feed back in next frame
%
%#codegen

LANDED_PHASE = uint8(5);
LATCHED = uint8(2);

% Allow reset only on ground with override active
if ground_override && (flight_state == LANDED_PHASE)
    current_response = uint8(0);
    active_response = current_response;
    return;
end

max_requested = uint8(0);

for i = 1:numel(fault_states)
    % If fault is latched, query its response tier from fault_config
    if fault_states(i).Persistence_State == LATCHED
        req = fault_config(i).Response_Tier;
    else
        req = uint8(0);
    end

    if req > max_requested
        max_requested = req;
    end
end

% Latch highest escalation level (prevent dropping back down mid-flight)
if max_requested > current_response
    current_response = max_requested;
end

active_response = current_response;
end
