"""Physics, configuration, and regression tests for the standalone ORV sizer."""

from __future__ import annotations

import json
import math
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np

ORV_DIRECTORY = Path(__file__).parents[1] / "isolated_sims" / "ORV_sizing"
sys.path.insert(0, str(ORV_DIRECTORY))

from ORV_config import DEFAULT_HOPPER_CONFIG, ORVConfigError, load_orv_config  # noqa: E402
from ORV_solver import (  # noqa: E402
    MM2_TO_M2,
    M2_TO_MM2,
    PSI_TO_PA,
    calculate_orv_relieving_temperature,
    calculate_sizing,
    compressible_mass_flow,
    solve_well_mixed_ullage_temperature,
)


class ORVSizingTests(unittest.TestCase):
    def test_choked_compressible_flow_matches_closed_form(self) -> None:
        upstream, temperature, cda = 1_000_000.0, 293.0, 2.5e-6
        gas_constant, gamma = 296.8, 1.4
        flow, choked = compressible_mass_flow(upstream, temperature, 100_000.0, cda, gas_constant, gamma)
        factor = math.sqrt(gamma / (gas_constant * temperature)) * (2.0 / (gamma + 1.0)) ** ((gamma + 1.0) / (2.0 * (gamma - 1.0)))
        self.assertTrue(bool(choked))
        self.assertAlmostEqual(float(flow), cda * upstream * factor, places=14)

    def test_unchoked_compressible_flow_is_selected_above_critical_ratio(self) -> None:
        flow, choked = compressible_mass_flow(1_000_000.0, 293.0, 900_000.0, 2.5e-6, 296.8, 1.4)
        self.assertFalse(bool(choked))
        self.assertGreater(float(flow), 0.0)

    def test_unit_conversion_constants_are_reciprocal_and_absolute(self) -> None:
        self.assertAlmostEqual(PSI_TO_PA, 6894.757293168)
        self.assertAlmostEqual(MM2_TO_M2 * M2_TO_MM2, 1.0)
        flow_psi, _ = compressible_mass_flow(1000.0 * PSI_TO_PA, 293.0, 14.8 * PSI_TO_PA, MM2_TO_M2, 296.8, 1.4)
        flow_pa, _ = compressible_mass_flow(6_894_757.293168, 293.0, 102_041.0, 1e-6, 296.8, 1.4)
        self.assertAlmostEqual(float(flow_psi), float(flow_pa), places=6)

    def test_pressure_sweep_and_linked_values(self) -> None:
        config = load_orv_config()
        result = calculate_sizing(config)
        self.assertGreaterEqual(len(result["pressure_psia"]), config["analysis"]["pressure_sweep_points"])
        self.assertTrue(np.all(np.diff(result["pressure_psia"]) > 0.0))
        self.assertAlmostEqual(result["maximum_inflow_pressure_psia"], config["tank"]["meop_psia"])
        self.assertAlmostEqual(result["maximum_gn2_inflow_kg_s"], config["failure"]["simulation_maximum_gn2_inflow_kg_s"])
        self.assertGreater(result["analytical_maximum_gn2_inflow_kg_s"], 0.0)
        self.assertAlmostEqual(result["gox"]["mass_flow_kg_s"], config["gox"]["provisional_generation_kg_s"])
        expected_flow = config["failure"]["gn2_mass_flow_safety_factor"] * result["maximum_gn2_inflow_kg_s"] + result["gox"]["mass_flow_kg_s"]
        self.assertAlmostEqual(result["combined_design_flow_kg_s"], expected_flow)
        self.assertAlmostEqual(result["governing_relieving_pressure_psia"], config["orv"]["maximum_relieving_pressure_psia"])
        self.assertGreater(result["required_gn2_cda_mm2"], 0.0)
        self.assertGreater(result["required_cda_mm2"], result["required_gn2_cda_mm2"])
        self.assertAlmostEqual(result["candidate_cda_mm2"], config["orv"]["candidate_cda_mm2"])

    def test_missing_required_hopper_input_is_explicit(self) -> None:
        raw = json.loads(DEFAULT_HOPPER_CONFIG.read_text(encoding="utf-8"))
        del raw["fluids"]["ambient_pressure_psia"]
        with tempfile.TemporaryDirectory() as directory:
            config_path = Path(directory) / "missing_ambient.json"
            config_path.write_text(json.dumps(raw), encoding="utf-8")
            with self.assertRaisesRegex(ORVConfigError, "ambient_pressure_psia"):
                load_orv_config(config_path)

    def test_ullage_energy_balance_conserves_energy(self) -> None:
        result = solve_well_mixed_ullage_temperature(0.25, 220.0, 0.03, 220.0, 90.0)
        self.assertAlmostEqual(result["temperature_k"], 220.0, places=7)
        self.assertAlmostEqual(result["net_heat_transfer_w"], 0.0, places=9)
        self.assertAlmostEqual(result["energy_residual_w"], 0.0, delta=1e-6)

    def test_ullage_temperature_inputs_are_kelvin(self) -> None:
        with self.assertRaisesRegex(ValueError, "kelvin"):
            solve_well_mixed_ullage_temperature(0.25, 0.0, 0.03, 150.0, 90.0)

    def test_limiting_thermal_cases(self) -> None:
        single_stream = solve_well_mixed_ullage_temperature(0.25, 293.0, 0.0, 150.0, 90.0)
        self.assertAlmostEqual(single_stream["temperature_k"], 293.0, places=6)
        strong_liquid_coupling = solve_well_mixed_ullage_temperature(
            0.25, 293.0, 0.03, 150.0, 90.0,
            liquid_to_ullage_ua_w_k=1.0e7,
        )
        self.assertAlmostEqual(strong_liquid_coupling["temperature_k"], 90.0, delta=0.02)

    def test_missing_thermal_inputs_retain_provisional_temperature(self) -> None:
        config = load_orv_config()
        gox = {"mass_flow_kg_s": 0.024847, "saturation_temperature_k": 147.0}
        thermal = calculate_orv_relieving_temperature(config, 0.293, gox)
        self.assertEqual(thermal["method"], "provisional_assumption_missing_thermal_closure")
        self.assertEqual(thermal["selected_temperature_k"], 150.0)
        self.assertIsNone(thermal["net_heat_transfer_w"])
        self.assertTrue(thermal["missing_inputs"])

    def test_energy_balance_numerical_convergence(self) -> None:
        result = solve_well_mixed_ullage_temperature(
            0.293, 293.0, 0.024847, 147.0, 90.0,
            liquid_to_ullage_ua_w_k=40.0,
            wall_temperature_k=250.0,
            wall_to_ullage_ua_w_k=20.0,
            net_external_sensible_heat_rate_w=100.0,
        )
        self.assertGreater(result["temperature_k"], 90.0)
        self.assertLess(result["temperature_k"], 293.0)
        self.assertAlmostEqual(result["energy_residual_w"], 0.0, delta=1e-5)


if __name__ == "__main__":
    unittest.main()
