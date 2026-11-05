"""Behavior checks for timing, numerical inputs, and run reporting."""

import csv
import json
import tempfile
import unittest
from pathlib import Path

import numpy as np
import h5py
from fullflow import Network, State, Lookup
from Hopper_Components import engine_efficiencies

from Hopper_Inputs import efficiency_value, table_value, ramp_profile, slew_position, delta_t_value
from Hopper_Control import SampledBangBang, ValveTravel, FunctionLookup
from Hopper_Report import sustained_maximum, summarize
from Hopper_Export import _safe_group_name
from Hopper_Elements import ClosedSafeDischargeCoefficient, Line


class InputTests(unittest.TestCase):
    def test_function_lookup_recomputes_as_pressure_changes(self):
        network = Network("test")
        pc = State(100 * 6894.757293168)
        lookup = FunctionLookup("efficiency", network, engine_efficiencies, chamber_pressure=pc,
                        mixture_ratio=1.8, outputs=("cstar", "cf"))
        self.assertAlmostEqual(lookup.cstar.value, .8)
        pc.value = 300 * 6894.757293168
        lookup.evaluate_states()
        self.assertAlmostEqual(lookup.cstar.value, .85)
        self.assertAlmostEqual(lookup.cf.value, .97)
    def test_tube_dimensions_are_consistent(self):
        line = Line("tube", outer_diameter=.01, wall_thickness=.001, length=1.)
        self.assertAlmostEqual(line.inner_diameter, .008)
        with self.assertRaises(ValueError):
            Line("bad tube", inner_diameter=.009, outer_diameter=.01, wall_thickness=.001, length=1.)
    def test_efficiency_csv_bilinear_and_endpoint_hold(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "eff.csv"
            path.write_text("pc_psia,mr,efficiency\n100,1,0.6\n300,1,0.8\n100,3,0.8\n300,3,1.0\n")
            self.assertAlmostEqual(table_value(path, 200*6894.757293168, 2), 0.8)
            self.assertAlmostEqual(table_value(path, 900*6894.757293168, 9), 1.0)
            self.assertAlmostEqual(efficiency_value(lambda p, mr: 0.3, path, 200*6894.757293168, 2), 0.8)

    def test_bad_efficiency_and_incomplete_grid_rejected(self):
        with self.assertRaises(ValueError):
            efficiency_value(lambda p, mr: 1.1, None, 1e6, 1.8)
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "bad.csv"
            path.write_text("pc_psia,mr,efficiency\n100,1,0.8\n300,1,0.8\n100,2,0.8\n")
            with self.assertRaises(ValueError):
                table_value(path, 1e6, 1.8)

    def test_ramp_endpoints_and_regen(self):
        profile = ramp_profile(400, 1600, 10, 101)
        self.assertEqual(profile[0], 400)
        self.assertEqual(profile[10], 1600)
        self.assertEqual(profile[5], 1000)
        self.assertEqual(delta_t_value([(0, 0), (2e6, 80)], 1e6), 40)


class TimingTests(unittest.TestCase):
    def test_closed_inertial_valve_does_not_divide_by_zero(self):
        network = Network("test")
        flow = State(.2)
        restriction = ClosedSafeDischargeCoefficient("closed", network, mass_flow=flow,
                upstream_pressure=2e6, downstream_pressure=1e6, density=800.,
                discharge_coefficient=0., cross_sectional_area=1e-5, length=.01)
        flow.store_previous()
        network.time.store_previous()
        network.time.value = .01
        restriction.evaluate_states()
        self.assertAlmostEqual(flow.previous + .01 * restriction.mass_flow_dot, 0.)
    def test_asymmetric_full_stroke_and_partial_reversal(self):
        self.assertAlmostEqual(slew_position(0, 1, 0.025, 0.05, 0.1), 0.5)
        self.assertAlmostEqual(slew_position(1, 0, 0.025, 0.05, 0.1), 0.75)
        self.assertAlmostEqual(slew_position(0.4, 0.1, 0.1, 0.05, 0.1), 0.1)

    def test_trial_evaluations_do_not_advance_valve_twice(self):
        network = Network("test")
        command, position = State(1.0), State(0.0)
        travel = ValveTravel("valve", network, command, position, 0.05, 0.10)
        position.store_previous()
        network.time.store_previous()
        network.time.value = 0.025
        travel.evaluate_states()
        self.assertAlmostEqual(position.value, 0.5)
        travel.evaluate_states()
        self.assertAlmostEqual(position.value, 0.5)
        # A rejected step retried at a smaller dt starts from the accepted state.
        network.time.value = 0.01
        travel.evaluate_states()
        self.assertAlmostEqual(position.value, 0.2)

    def test_controller_holds_between_samples_and_preserves_hysteresis(self):
        network = Network("test")
        pressure, command = State(100.0), State(0.0)
        control = SampledBangBang("controller", network, pressure, command, 100., 5., 5., 0.05)
        for state in (network.time, pressure, command, control.sample_index):
            state.store_previous()
        control.evaluate_states()
        command.store_previous()
        control.sample_index.store_previous()
        pressure.value = 90
        pressure.store_previous()
        network.time.value = 0.025
        network.time.store_previous()
        control.evaluate_states()
        self.assertEqual(command.value, 0)
        network.time.value = 0.05
        network.time.store_previous()
        control.evaluate_states()
        self.assertEqual(command.value, 1)


class ReportTests(unittest.TestCase):
    def test_saved_run_reports_capacity_residuals_and_copv_limits(self):
        time = np.arange(0, 2.01, .1)
        tracks = {"Thrust [N]": np.full_like(time, 1000.),
                  "Target Thrust [N]": np.full_like(time, 1100.),
                  "Chamber Pressure [psia]": np.full_like(time, 200.),
                  "Ox Mass Flow [kg/s]": np.full_like(time, .9),
                  "Fuel Mass Flow [kg/s]": np.full_like(time, .5),
                  "LOX Liquid Mass [kg]": 10 - time,
                  "Fuel Liquid Mass [kg]": 5 - time / 2,
                  "OTV CdA [mm2]": np.full_like(time, 20.),
                  "COPV 1 Temperature [K]": 293 - 3*time}
        metadata = {"profile": "test", "profile_start": 0., "profile_end": 2., "hold_time": .5,
                    "tracking_tolerance": .02, "map_pc_max_pa": 1e7, "map_mr_min": .5, "map_mr_max": 4.,
                    "tanks": {"LOX": {"initial_mass": 10.}, "Fuel": {"initial_mass": 5.}},
                    "throttles": {"LOX": [{"name": "OTV", "max_cda": 20e-6}]},
                    "copvs": [{"name": "COPV 1", "minimum_temperature": 290.}]}
        with tempfile.TemporaryDirectory() as directory:
            filename = Path(directory) / "report.h5"
            with h5py.File(filename, "w") as h5:
                h5.attrs["hopper_metadata"] = json.dumps(metadata)
                run = h5.create_group("test/transient/runs/base")
                run.create_dataset("time", data=time)
                saved = run.create_group("tracks")
                for name, values in tracks.items():
                    saved.create_dataset(_safe_group_name(name), data=values)
            report = summarize(filename)
            self.assertEqual(report["demonstrated_maximum_N"], 1000.)
            self.assertTrue(report["capacity_shortfall_demonstrated"])
            self.assertEqual(report["residuals"]["LOX"]["final_kg"], 8.)
            self.assertEqual(report["residuals"]["LOX"]["remaining_percent"], 80.)
            self.assertTrue(report["copv_temperatures"][0]["violated"])

    def test_spike_is_not_sustained_thrust(self):
        time = np.arange(0, 2.01, 0.1)
        thrust = np.full_like(time, 1000.)
        thrust[10] = 2000
        value, interval = sustained_maximum(time, thrust, np.ones(len(time), bool), 0.5)
        self.assertEqual(value, 1000)
        self.assertGreaterEqual(interval[1] - interval[0], 0.5 - 1e-9)

    def test_invalid_samples_break_sustained_interval(self):
        time = np.arange(0, 0.81, 0.1)
        valid = np.ones(len(time), bool)
        valid[4] = False
        self.assertEqual(sustained_maximum(time, time, valid, 0.5), (None, None))


if __name__ == "__main__":
    unittest.main()
