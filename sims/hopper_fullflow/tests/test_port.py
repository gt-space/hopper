"""Checks for the pieces brought over from the element-list model: throttle slew limit,
ramp breakpoints, relief opening, and the bottle / header energy balance."""
import tempfile
import unittest
from pathlib import Path

import h5py
import numpy as np
from fullflow import CompressibleOrifice, Network, State, Transient, Volume

import Backend.Hopper_Prop_System as simulation
from Backend.Hopper_Elements import ReliefValve


class ScheduleTests(unittest.TestCase):
    def test_step_duration_follows_the_slew_limit_within_its_bounds(self):
        limit, step = simulation.THROTTLE_RATE_LIMIT, simulation.THROTTLE_STEP_TIME
        duration = simulation.step_duration({"time": 0., "thrust": 1000.}, {"time": 5., "thrust": 1000. + limit})
        self.assertAlmostEqual(duration, 1.)
        small = simulation.step_duration({"time": 0., "thrust": 1000.}, {"time": 5., "thrust": 1000.})
        self.assertAlmostEqual(small, step)
        short_gap = simulation.step_duration({"time": 0., "thrust": 1000.}, {"time": .1, "thrust": 1000. + limit})
        self.assertAlmostEqual(short_gap, .09)

    def test_ramps_have_breakpoints_and_end_closed(self):
        targets = [{"time": 0., "thrust": 1000., "ox_mass_flow": .5}, {"time": 1., "thrust": 1000., "ox_mass_flow": .5}]
        times, values = simulation.target_schedule(targets, "ox_mass_flow", "step")
        self.assertTrue(all(b > a for a, b in zip(times, times[1:])))
        ramp = simulation.PROFILE_START_TIME - simulation.MAIN_OX_OPEN_TIME
        opening = [t for t in times if simulation.MAIN_OX_OPEN_TIME < t <= simulation.PROFILE_START_TIME]
        self.assertGreaterEqual(len(opening), ramp / simulation.RAMP_BREAKPOINT_SPACING - 1)
        self.assertEqual(values[0], 0.)
        self.assertAlmostEqual(values[-1], 0.)


class ReliefTests(unittest.TestCase):
    def test_relief_opens_linearly_from_crack_to_full(self):
        valve = ReliefValve("relief", 1e-6, crack_pressure=100.)
        self.assertAlmostEqual(valve.full_open_pressure, 110.)
        upstream = State(100.)
        opening = valve.relief_opening(upstream, 0.)
        self.assertAlmostEqual(float(opening.value), 0.)
        upstream.value = 105.
        self.assertAlmostEqual(float(opening.value), .5)
        upstream.value = 500.
        self.assertAlmostEqual(float(opening.value), 1.)


class HeaderTests(unittest.TestCase):
    def test_header_energy_uses_enthalpy_per_kilogram(self):
        """Gas leaving the bottles should reach the header near bottle temperature. Passing the
        enthalpy flow as total_enthalpy_in scales the incoming energy by the flow rate and
        the header temperature runs away from the bottles'."""
        saved = simulation.COPV_MODEL
        simulation.COPV_MODEL = "header"
        try:
            network = Network("Header Check")
            header, gas, _ = simulation.build_copvs(network)
            vent = CompressibleOrifice("Vent", network, upstream_total_pressure=header.pressure,
                                       upstream_total_temperature=gas.temperature, downstream_pressure=simulation.AMBIENT_PRESSURE,
                                       discharge_coefficient=1.0, cross_sectional_area=1e-6, gas_constant=gas.gas_constant,
                                       specific_heat_ratio=gas.gamma)
            header.mass_flow_out = State.maximum(vent.mass_flow, 0.0)
            network.track("Header T", gas.temperature)
            with tempfile.TemporaryDirectory() as directory:
                filename = str(Path(directory) / "header.h5")
                Transient(network).solve(dt=.01, t_final=.1, save_dt=.01, filename=filename, rtol=1e-6,
                                         verbose=False, statistics=False)
                with h5py.File(filename, "r") as h5:
                    tracks = h5["Header_Check/transient/runs/base/tracks"]
                    header_t = np.asarray(tracks["Header_T"][()], float)
                    bottle_t = np.asarray(tracks[f"{simulation.COPVS[0].name.replace(' ', '_')}_Temperature_[K]"][()], float)
            self.assertLess(abs(header_t[-1] - bottle_t[-1]), 5.0)
        finally:
            simulation.COPV_MODEL = saved


if __name__ == "__main__":
    unittest.main()
