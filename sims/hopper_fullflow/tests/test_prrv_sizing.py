"""Validation tests for the standalone high-pressure PRRV sizing screen."""

from __future__ import annotations

import json
import copy
import math
import sys
import tempfile
import unittest
from pathlib import Path

PRRV_DIRECTORY = Path(__file__).parents[1] / "isolated_sims" / "PRRV_sizing"
sys.path.insert(0, str(PRRV_DIRECTORY))

from PRRV_config import DEFAULT_HOPPER_CONFIG, PRRVConfigError, load_prrv_config  # noqa: E402
from PRRV_solver import MM2_TO_M2, calculate, pressure_to_psia  # noqa: E402


class PRRVSizingTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.config = load_prrv_config()
        cls.result = calculate(cls.config)

    def test_repository_candidate_and_fill_load_status_are_reused(self) -> None:
        self.assertIsNone(self.config["failure"]["maximum_fill_mass_flow_kg_s"])
        self.assertEqual(self.config["candidate"]["configured_cda_mm2"], 5.0)
        self.assertEqual(self.config["candidate"]["configured_crack_pressure_psig"], 4950.0)

    def test_pressure_reference_conversion(self) -> None:
        self.assertEqual(pressure_to_psia(4550.0, "psia", 14.8), 4550.0)
        self.assertEqual(pressure_to_psia(4550.0, "psig", 14.8), 4564.8)
        with self.assertRaisesRegex(ValueError, "resolved"):
            pressure_to_psia(4550.0, "unconfirmed", 14.8)

    def test_mass_flow_factor_is_applied_once_in_trade_mode(self) -> None:
        for inflow, design in zip(self.result["trial_fill_flow_kg_s"], self.result["trial_design_flow_kg_s"]):
            self.assertAlmostEqual(design, self.config["failure"]["mass_flow_factor"] * inflow, places=12)

    def test_real_gas_required_area_reproduces_each_trade_flow(self) -> None:
        reproduced = self.result["trial_required_cda_mm2"] * MM2_TO_M2 * self.result["real_gas"]["mass_flux_kg_m2_s"]
        self.assertTrue(all(math.isclose(a, b) for a, b in zip(reproduced, self.result["trial_design_flow_kg_s"])))
        self.assertTrue(self.result["real_gas"]["is_choked"])

    def test_high_pressure_real_gas_effect_is_reported(self) -> None:
        self.assertGreater(self.result["real_gas"]["upstream_compressibility_factor"], 1.0)
        self.assertNotAlmostEqual(self.result["real_gas"]["mass_flux_kg_m2_s"], self.result["ideal_mass_flux_kg_m2_s"], places=3)

    def test_configured_setpoint_and_candidate_status(self) -> None:
        expected_setpoint = self.config["candidate"]["configured_crack_pressure_psig"] + self.config["ambient_pressure_psia"]
        self.assertEqual(self.config["pressure"]["set_value_psi"], expected_setpoint)
        self.assertEqual(self.config["pressure"]["set_reference"], "psia")
        self.assertEqual(self.result["evaluation_pressure_psia"], expected_setpoint)
        self.assertEqual(self.result["candidate_analytical_status"], "UNRESOLVED_MAXIMUM_FLOW")
        self.assertIn("UNVERIFIED", self.result["manufacturer_capacity_validation_status"])

    def test_known_flow_mode_is_supported(self) -> None:
        known = copy.deepcopy(self.config)
        known["failure"]["maximum_fill_mass_flow_kg_s"] = 0.2
        result = calculate(known)
        self.assertEqual(result["mode"], "A_known_maximum_fill_flow")
        self.assertAlmostEqual(result["design_flow_kg_s"], 0.3)
        self.assertAlmostEqual(result["required_cda_mm2"] * MM2_TO_M2 * result["real_gas"]["mass_flux_kg_m2_s"], 0.3)

    def test_missing_candidate_data_is_explicit(self) -> None:
        raw = json.loads(DEFAULT_HOPPER_CONFIG.read_text(encoding="utf-8"))
        prrv = next(branch for branch in raw["side_paths"] if branch["name"] == "PRRV Relief")
        del prrv["elements"][0]["cda_mm2"]
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "missing_candidate.json"
            path.write_text(json.dumps(raw), encoding="utf-8")
            with self.assertRaisesRegex(PRRVConfigError, "cda_mm2"):
                load_prrv_config(path)


if __name__ == "__main__":
    unittest.main()
