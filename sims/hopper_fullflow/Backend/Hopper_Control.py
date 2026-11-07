"""Retry-safe valve travel and sampled hysteresis using accepted FullFlow states."""

import math

from fullflow import Component, State, Lookup
from Backend.Hopper_Inputs import slew_position


def FunctionLookup(name, network, function, **inputs):
    """A function must be called again, rather than mutating its previous dict result.

    FullFlow 3 treats dict.update as an object updater by default. Disabling
    object reuse keeps Pc-dependent calculated outputs live. This is an attribute
    setting in FullFlow 3, not a constructor keyword (keywords go to function).
    """
    lookup = Lookup(name, network, function, **inputs)
    lookup.reuse_existing = False
    return lookup


def previous(state):
    try:
        return state.previous
    except ValueError:
        return state.value


class ValveTravel(Component):
    """Opening/closing times are linear full-stroke travel times in seconds."""

    def __init__(self, name, network, command, position, open_time=0.0, close_time=0.0):
        if min(open_time, close_time) < 0 or not math.isfinite(open_time + close_time):
            raise ValueError(f"{name}: valve travel times must be finite and nonnegative.")
        self.setup()

    def evaluate_states(self):
        dt = max(0.0, float(self.network.time.value) - float(previous(self.network.time)))
        self.position.value = slew_position(previous(self.position), self.command.value, dt,
                                            self.open_time.value, self.close_time.value)


class SampledBangBang(Component):
    """Evaluate pressure at sample instants; hold the last command between samples.

    Reads accepted states, never solver trial pressures. Sampling is recognized at
    the first accepted step at/after a sample instant. Integration step is capped
    by the caller; remaining timing quantization is at most one integration step.
    """

    def __init__(self, name, network, pressure, command, set_pressure, upper_limit,
                 lower_limit, sample_period, horizon=60.0):
        if (not all(math.isfinite(x) for x in (sample_period, set_pressure, upper_limit, lower_limit, horizon))
                or sample_period <= 0 or min(upper_limit, lower_limit) < 0 or lower_limit >= set_pressure):
            raise ValueError(f"{name}: sample period must be positive; deadbands nonnegative.")
        self.setup()
        self.sample_index = State(-1.0)

    def command_breakpoints(self):
        return tuple(i * self.sample_period.value for i in
                     range(1, math.ceil(self.horizon.value / self.sample_period.value) + 1))

    def evaluate_states(self):
        start_time = float(previous(self.network.time))
        index = math.floor((start_time + 1e-10) / self.sample_period.value)
        command = float(previous(self.command))
        self.sample_index.value = float(previous(self.sample_index))
        if index > previous(self.sample_index):
            p = float(previous(self.pressure))
            if p > self.set_pressure.value + self.upper_limit.value:
                command = 0.0
            elif p < self.set_pressure.value - self.lower_limit.value:
                command = 1.0
            self.sample_index.value = float(index)
        self.command.value = command
