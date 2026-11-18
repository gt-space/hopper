"""All relief-sizer inputs must follow the shared Hopper JSON."""

from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

from isolated_sims.FRV_sizing.frv_config import load_frv_config
from isolated_sims.ORV_sizing.ORV_config import load_orv_config
from isolated_sims.PRRV_sizing.PRRV_config import load_prrv_config


ROOT = Path(__file__).parents[1]
DEFAULT_CONFIG = ROOT / "hopper_config.json"


class ReliefSizingConfigLinkTests(unittest.TestCase):
    def _write(self, raw: dict) -> tuple[tempfile.TemporaryDirectory, Path]:
        directory = tempfile.TemporaryDirectory()
        path = Path(directory.name) / "hopper_config.json"
        path.write_text(json.dumps(raw), encoding="utf-8")
        return directory, path

    def test_current_and_new_pid_names_are_supported_and_values_stay_linked(self) -> None:
        raw = json.loads(DEFAULT_CONFIG.read_text(encoding="utf-8"))
        raw["relief_sizing"]["shared"]["mass_flow_safety_factor"] = 1.6
        raw["relief_sizing"]["shared"]["pressure_sweep_points"] = 123
        raw["tanks"]["fuel"]["set_pressure_psia"] = 510.0
        raw["tanks"]["oxidizer"]["set_pressure_psia"] = 515.0

        fuel_path = raw["flow_paths"]["fuel_pressurant"]
        fuel_orifice = next(item for item in fuel_path if item["name"] == "OR-2")
        fuel_orifice["name"] = "FBBORF"
        fuel_orifice["cda_mm2"] = 4.25
        next(item for item in fuel_path if item["name"] == "PRFCK")["cda_mm2"] = 44.2
        oxidizer_path = raw["flow_paths"]["oxidizer_pressurant"]
        oxidizer_orifice = next(item for item in oxidizer_path if item["name"] == "OR-1")
        oxidizer_orifice["name"] = "OBBORF"
        oxidizer_orifice["cda_mm2"] = 8.75
        next(item for item in oxidizer_path if item["name"] == "PROCK")["cda_mm2"] = 44.2
        raw["regulators"][0]["name"] = "PUREG"

        frv = next(branch for branch in raw["side_paths"] if branch["name"] == "FRV Relief")["elements"][0]
        orv = next(branch for branch in raw["side_paths"] if branch["name"] == "ORV Relief")["elements"][0]
        frv["crack_pressure_psig"] = 585.2
        frv["cda_mm2"] = 16.0
        orv["crack_pressure_psig"] = 585.2
        orv["cda_mm2"] = 31.0

        directory, path = self._write(raw)
        self.addCleanup(directory.cleanup)
        fuel = load_frv_config(path)
        oxidizer = load_orv_config(path)

        self.assertEqual(fuel["failure"]["orifice_name"], "FBBORF")
        self.assertEqual(fuel["failure"]["orifice_cda_mm2"], 4.25)
        self.assertEqual(fuel["failure"]["check_valve_cda_mm2"], 44.2)
        self.assertEqual(fuel["context"]["regulator_name"], "PUREG")
        self.assertEqual(fuel["tank"]["meop_psia"], 510.0)
        self.assertEqual(fuel["frv"]["proposed_set_pressure_psia"], 600.0)
        self.assertEqual(fuel["frv"]["configured_full_open_cda_mm2"], 16.0)
        self.assertAlmostEqual(fuel["failure"]["sizing_margin_fraction"], 0.6)
        self.assertEqual(fuel["analysis"]["pressure_sweep_points"], 123)

        self.assertEqual(oxidizer["failure"]["failing_orifice_name"], "OBBORF")
        self.assertEqual(oxidizer["failure"]["orifice_cda_mm2"], 8.75)
        self.assertEqual(oxidizer["failure"]["check_valve_cda_mm2"], 44.2)
        self.assertEqual(oxidizer["tank"]["meop_psia"], 515.0)
        self.assertEqual(oxidizer["orv"]["proposed_set_pressure_psia"], 600.0)
        self.assertEqual(oxidizer["orv"]["candidate_cda_mm2"], 31.0)
        self.assertEqual(oxidizer["failure"]["gn2_mass_flow_safety_factor"], 1.6)
        self.assertEqual(oxidizer["analysis"]["pressure_sweep_points"], 123)

    def test_prrv_pressure_and_temperature_follow_hardware(self) -> None:
        raw = json.loads(DEFAULT_CONFIG.read_text(encoding="utf-8"))
        raw["bottles"][0]["pressure_psia"] = 4600.0
        raw["bottles"][0]["temperature_k"] = 285.0
        raw["relief_sizing"]["shared"]["mass_flow_safety_factor"] = 1.7
        prrv = next(branch for branch in raw["side_paths"] if branch["name"] == "PRRV Relief")["elements"][0]
        prrv["crack_pressure_psig"] = 4985.2

        directory, path = self._write(raw)
        self.addCleanup(directory.cleanup)
        config = load_prrv_config(path)

        self.assertEqual(config["pressure"]["meop_value_psi"], 4600.0)
        self.assertEqual(config["pressure"]["set_value_psi"], 5000.0)
        self.assertEqual(config["relief"]["temperature_k"], 285.0)
        self.assertEqual(config["failure"]["mass_flow_factor"], 1.7)


if __name__ == "__main__":
    unittest.main()
