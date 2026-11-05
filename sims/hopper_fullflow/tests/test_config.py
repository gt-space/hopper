"""Config migration, validation and reproducible run archives."""
import copy
import json
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest

from Hopper_Config import load_config, build_component, efficiency, integration_step, ROOT
from Hopper_Run_Config import prepare_run


class ConfigTests(unittest.TestCase):
    def test_default_geometry_and_timing_preserved(self):
        config, data, _ = load_config(ROOT/'hopper_config.json')
        self.assertAlmostEqual(data['THROAT_AREA'], .00067226)
        self.assertAlmostEqual(data['LOX_FEED'][0].inner_diameter, (.375-2*.035)*.0254)
        self.assertAlmostEqual(integration_step(data), .005)
        self.assertEqual(len(data['COPVS']), len(config['bottles']))
        self.assertEqual(data['LOX_TANK'].ullage_htc, config['tanks']['oxidizer']['ullage_htc_w_m2_k'])

    def test_unknown_and_duplicate_fields_rejected(self):
        config, _, _ = load_config(ROOT/'hopper_config.json')
        with tempfile.TemporaryDirectory() as directory:
            file = Path(directory)/'config.json'
            config['simulation']['maximun_step_s'] = .02
            file.write_text(json.dumps(config))
            with self.assertRaisesRegex(ValueError,'unknown fields'):
                load_config(file)
            file.write_text('{"schema_version":1,"schema_version":1}')
            with self.assertRaisesRegex(ValueError,'Duplicate'):
                load_config(file)

    def test_legacy_config_gets_bounded_solver_default(self):
        config, _, _ = load_config(ROOT/'hopper_config.json')
        config['simulation'].pop('nonlinear_max_evaluations')
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/'config.json'
            path.write_text(json.dumps(config))
            resolved, data, _ = load_config(path)
            self.assertEqual(data['PROP_MAX_NFEV'], 20)
            self.assertEqual(resolved['simulation']['nonlinear_max_evaluations'], 20)
            for invalid in (0, -1, 1.5, True):
                config['simulation']['nonlinear_max_evaluations'] = invalid
                path.write_text(json.dumps(config))
                with self.assertRaisesRegex(ValueError, 'nonlinear_max_evaluations'):
                    load_config(path)

    def test_geometry_and_negative_valve_times_rejected(self):
        with self.assertRaises(ValueError):
            build_component({'type':'Line','name':'bad','outer_diameter_in':.375,
                'wall_thickness_in':.035,'inner_diameter_in':.1,'length_in':12})
        with self.assertRaises(ValueError):
            build_component({'type':'Valve','name':'bad','cda_mm2':2,'open_time_s':-.1})

    def test_defaults_are_materialized(self):
        spec={'type':'Valve','name':'test','cda_mm2':2}
        build_component(spec)
        self.assertEqual(spec['close_time_s'],0)
        self.assertIn('schedule',spec)

    def test_function_and_relative_csv(self):
        function,_ = efficiency({'mode':'polynomial','coefficients':[[.8,.01],[.0001]]}, ROOT)
        self.assertAlmostEqual(function(100*6894.757293168361,2),.83)
        with tempfile.TemporaryDirectory() as directory:
            file=Path(directory)/'eff.csv'
            file.write_text('pc_psia,efficiency\n100,0.8\n300,0.9\n')
            _,path=efficiency({'mode':'csv','csv':'eff.csv'},Path(directory))
            self.assertEqual(Path(path),file)

    def test_run_archive_remains_replayable_and_independent(self):
        config, data, path=load_config(ROOT/'hopper_config.json')
        with tempfile.TemporaryDirectory() as directory:
            source=Path(directory)/'eff.csv'
            source.write_text('pc_psia,efficiency\n100,0.8\n300,0.9\n')
            config['engine']['cstar_efficiency']={'mode':'csv','csv':str(source)}
            data['CSTAR_EFFICIENCY_CSV']=str(source)
            data['PROP_RESULTS_DIRECTORY']=directory
            namespace=SimpleNamespace(**data,CONFIG=config,CONFIG_PATH=path,CONFIG_SOURCE_TEXT=json.dumps(config))
            first=prepare_run(namespace)
            source.write_text('changed after run')
            archived, replay, _=load_config(first/'run_config.json')
            self.assertEqual(archived['engine']['cstar_efficiency']['csv'],'inputs/cstar_efficiency.csv')
            self.assertTrue(Path(replay['CHAMBER_MAP_FILENAME']).is_file())
            self.assertIn('0.8',Path(replay['CSTAR_EFFICIENCY_CSV']).read_text())
            self.assertAlmostEqual(json.loads((first/'resolved_inputs.json').read_text())['derived']['integration_step_cap_s'],.005)
            self.assertTrue((first/'input_config.json').exists())
            namespace.PROP_RESULTS_DIRECTORY=directory
            second=prepare_run(namespace)
            self.assertNotEqual(first,second)
            self.assertTrue((first/'run_config.json').exists())


if __name__ == '__main__':
    unittest.main()
