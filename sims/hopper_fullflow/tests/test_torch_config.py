"""Topology, shared purge tee, unset measurements and torch-relative OTHX inputs."""
import copy
import json
import tempfile
import unittest
from pathlib import Path

from Backend.Hopper_Config import load_config, ROOT
from Backend.Hopper_Torch_Config import heat_exchanger_outlet_temperature, validate_time_curve


class TorchConfigTests(unittest.TestCase):
    def setUp(self):
        self.config = json.loads((ROOT/'hopper_config.json').read_text())
        # Temporary configs retain valid external asset paths.
        for entry in self.config['thrust_profile']['profiles'].values():
            if 'csv' in entry:
                entry['csv'] = str(ROOT/entry['csv'])

    def load(self, config=None):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/'torch.json'
            path.write_text(json.dumps(config or self.config))
            return load_config(path)

    def test_requested_order_shared_chamber_and_purge_tee(self):
        _, data, _ = self.load()
        net = data['TORCH_NETWORK']
        fuel = net['fuel_leg']
        ox = net['oxidizer_leg']
        purge = net['purge_leg']
        self.assertEqual(fuel['source'], 'FTV Inlet')
        self.assertEqual([e['name'] for e in fuel['elements']], ['FIMV','FICK','FIORF','Torch Fuel Injector'])
        self.assertEqual(ox['source'], 'OTV Inlet')
        self.assertEqual([e['name'] for e in ox['elements']], ['OIMV','OICK','OIORF','OTHX','Torch Ox Injector'])
        self.assertEqual([e['name'] for e in purge['elements']], ['TFPV','FIPCK','FIPORF'])
        main_purge = next(b for b in data['BRANCHES'] if b.name == 'Fuel Purge')
        self.assertEqual(purge['source'], main_purge.source)
        self.assertEqual(purge['sink'], fuel['elements'][2]['outlet_node'])
        self.assertEqual(purge['sink'], 'FIORF Outlet')
        self.assertEqual(fuel['sink'], net['chamber']['name'])
        self.assertEqual(ox['sink'], net['chamber']['name'])
        self.assertEqual(net['throat']['source'], net['chamber']['name'])
        self.assertEqual(net['throat']['sink'], 'Main Chamber')
        self.assertFalse(data['TORCH_ENABLED'])
        self.assertFalse(any(b.name == 'Torch Fuel Purge' for b in data['BRANCHES']))
        self.assertEqual(data['TORCH_BRANCHES'], [])

    def test_unknown_hardware_is_explicitly_unset(self):
        # Unset fields must remain loadable for configuration work, independent
        # of the user's subsequently populated default hardware measurements.
        net = self.config['torch']['network']
        net['chamber']['volume_l'] = None
        net['throat']['area_in2'] = None
        net['oxidizer_leg']['elements'][3]['length_in'] = None
        net['fuel_leg']['elements'][2]['cda_mm2'] = None
        _, data, _ = self.load()
        net = data['TORCH_NETWORK']
        self.assertIsNone(net['chamber']['volume_l'])
        self.assertIsNone(net['throat']['area_in2'])
        self.assertIsNone(net['oxidizer_leg']['elements'][3]['length_in'])
        self.assertIsNone(net['fuel_leg']['elements'][2]['cda_mm2'])

    def test_hardware_inventory_includes_configured_torch_parts_as_unsolved(self):
        from Backend.Hopper_Hardware import hardware_rows
        rows = {row[1]: row for row in hardware_rows()}
        for name in ('FIMV','FICK','FIORF','TFPV','FIPCK','FIPORF','OIMV','OICK','OIORF','OTHX',
                     'Torch Fuel Injector','Torch Ox Injector','Torch Chamber','Torch Throat'):
            self.assertIn(name, rows)
            self.assertIn('transient coupling pending', rows[name][4])

    def test_othx_temperature_rise_uses_torch_open_time(self):
        hx = copy.deepcopy(self.config['torch']['network']['oxidizer_leg']['elements'][3])
        hx['temperature_rise']['time_s_delta_t_k'] = [[0.,10.],[.2,50.]]
        self.assertEqual(heat_exchanger_outlet_temperature(hx,90.,.19,.2),90.)
        self.assertEqual(heat_exchanger_outlet_temperature(hx,90.,.2,.2),100.)
        self.assertAlmostEqual(heat_exchanger_outlet_temperature(hx,90.,.3,.2),120.)
        self.assertEqual(heat_exchanger_outlet_temperature(hx,90.,1.,.2),140.)

    def test_bad_time_curves_rejected(self):
        for points in ([], [[0,0],[0,1]], [[.2,1],[.1,2]], [[-1,0]], [[0,float('nan')]], [[0]], [[True,0]]):
            with self.subTest(points=points), self.assertRaises(ValueError):
                validate_time_curve(points)

    def test_bad_connectivity_and_geometry_rejected(self):
        mutations = (
            lambda n: n['purge_leg'].update(sink='not a fuel node'),
            lambda n: n['oxidizer_leg'].update(source='not a supply'),
            lambda n: n['fuel_leg'].update(sink='Main Chamber'),
            lambda n: n['throat'].update(source='FIORF Outlet'),
            lambda n: n['throat'].update(area_in2=-1),
            lambda n: n['oxidizer_leg']['elements'][3].update(outer_diameter_in=.25,wall_thickness_in=.2),
            lambda n: n['fuel_leg']['elements'][0].update(cda_mm2=-1),
            lambda n: n['purge_leg']['elements'][0].update(name='FIMV'),
        )
        for mutation in mutations:
            config = copy.deepcopy(self.config)
            mutation(config['torch']['network'])
            with self.assertRaises(ValueError):
                self.load(config)

    def test_enabling_unimplemented_network_fails_instead_of_using_legacy_shortcut(self):
        self.config['torch']['enabled'] = True
        with self.assertRaisesRegex(ValueError,'transient builder is not implemented'):
            self.load()

    def test_previous_schema_two_torch_branches_still_load(self):
        self.config['torch'].pop('network')
        self.config['torch']['branches'] = [
            {'name':'Legacy Torch Fuel','source':'FTV Inlet','sink':'Main Chamber','into_sink':True,
             'elements':[{'type':'Valve','name':'Legacy FIMV','cda_mm2':2,'schedule':'torch'}]}]
        _,data,_ = self.load()
        self.assertIsNone(data['TORCH_NETWORK'])
        self.assertEqual(data['TORCH_BRANCHES'][0].name,'Legacy Torch Fuel')


if __name__ == '__main__':
    unittest.main()
