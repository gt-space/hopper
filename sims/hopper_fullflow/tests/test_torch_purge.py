"""Independent torch purge physics, event scheduling, backflow and conservation."""
import copy
import unittest

import numpy as np
from scipy.constants import psi
from fullflow import CompressibleOrifice, Network

from Backend.Hopper_Config import load_config, ROOT
from Backend.Hopper_Torch_Purge import TorchPurgeModel, gas_chain, orifice_flow, valve_position
from Backend.Hopper_Torch_Config import resolve_purge_sequence
from Backend.Hopper_Torch_Diagram import torch_connections


class GasPathTests(unittest.TestCase):
    def test_orifice_matches_fullflow_for_choked_and_subcritical_flow(self):
        for back in (1e5, 4.5e5):
            network=Network('reference')
            part=CompressibleOrifice('orifice',network,upstream_total_pressure=5e5,
                downstream_pressure=back,upstream_total_temperature=293.,cross_sectional_area=2e-6,
                discharge_coefficient=1.,gas_constant=296.8,specific_heat_ratio=1.4)
            part.evaluate_states()
            self.assertAlmostEqual(orifice_flow(2e-6,5e5,back,293.,1.4,296.8),part.mass_flow.value,places=11)

    def test_each_series_restriction_has_the_same_flow(self):
        areas=[2.2e-6,5e-6,8e-6,.1051e-6]
        flow,nodes=gas_chain(areas,764.8*psi,300.*psi,293.,1.4,296.8)
        self.assertGreater(flow,0)
        self.assertTrue(np.all(np.diff(nodes)<0))
        for a,up,down in zip(areas,nodes,nodes[1:]):
            self.assertAlmostEqual(orifice_flow(a,up,down,293.,1.4,296.8),flow,places=10)
        self.assertAlmostEqual(nodes[0],764.8*psi)
        self.assertAlmostEqual(nodes[-1],300.*psi)

    def test_closed_and_reverse_supply_flow_are_blocked(self):
        flow,nodes=gas_chain([0.,1e-6],5e5,1e5,293.,1.4,296.8)
        self.assertEqual(flow,0)
        self.assertEqual(nodes,[5e5,1e5,1e5])
        self.assertEqual(gas_chain([1e-6],1e5,5e5,293.,1.4,296.8)[0],0.)

    def test_very_small_opening_and_choked_endpoint_roundoff(self):
        for area in (1e-10,1e-9,1e-8,1e-6):
            flow,nodes=gas_chain([area,5e-6,.1051e-6],764.8*psi,150.*psi,293.,1.4,296.8)
            self.assertGreater(flow,0)
            self.assertTrue(np.all(np.isfinite(nodes)))


class ValveCommandTests(unittest.TestCase):
    def test_travel_and_partial_reversal(self):
        valve={'command_events_s':[[0.,0.],[1.,1.],[1.04,0.]],'open_time_s':.08,'close_time_s':.1}
        self.assertEqual(valve_position(valve,1.),0.)
        self.assertAlmostEqual(valve_position(valve,1.02),.25)
        self.assertAlmostEqual(valve_position(valve,1.04),.5)
        self.assertAlmostEqual(valve_position(valve,1.065),.25)
        self.assertEqual(valve_position(valve,1.2),0.)

    def test_zero_travel_is_instant_at_command_event(self):
        valve={'command_events_s':[[1.,1.],[2.,0.]],'open_time_s':0.,'close_time_s':0.}
        self.assertEqual(valve_position(valve,.9),0.)
        self.assertEqual(valve_position(valve,1.),1.)
        self.assertEqual(valve_position(valve,2.),0.)


class PurgeTransientTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.config,cls.values,_=load_config(ROOT/'hopper_config.json')
        cls.config=copy.deepcopy(cls.config)
        # Short legacy-format case also verifies that saved per-valve archives replay.
        settings,network=resolve_purge_sequence(cls.config['torch']['standalone_purge'],cls.config['torch']['network'])
        cls.config['torch']['standalone_purge']=settings
        cls.config['torch']['network']=network
        settings=cls.config['torch']['standalone_purge']
        # Fixed cold-gas regression hardware: editable vehicle inputs must not
        # change this test's reference pressure margin or flow expectations.
        net=cls.config['torch']['network']
        net['chamber'].update(volume_l=.0019,initial_pressure_psia=150.,initial_temperature_k=293.)
        net['throat']['area_in2']=.006504
        for part,area in zip(net['purge_leg']['elements'],(2.2,5.,8.)):
            part['cda_mm2']=area
        net['fuel_leg']['elements'][-1]['cda_mm2']=.1051
        settings.update(purge_supply_pressure_psia=764.8,purge_supply_temperature_k=293.,
                        main_chamber_gas_temperature_k=293.,gamma=1.4)
        settings.update(imposed_pressure_until_s=.02,final_time_s=.25,maximum_step_s=.005,
                        save_interval_s=.002,ladder_times_s=[.01,.02,.05,.1,.2],
                        main_chamber_pressure_time_s_psia=[[0.,14.8],[.02,300.]])
        for key in ('fuel_leg','oxidizer_leg','purge_leg'):
            for part in cls.config['torch']['network'][key]['elements']:
                if part['type']=='Valve':
                    part['command_events_s']=([[0.,0.],[.02,1.],[.18,0.]] if key=='purge_leg'
                                              else [[0.,1.],[.02,0.]])
                    part['open_time_s']=.008
                    part['close_time_s']=.008
        cls.model=TorchPurgeModel(cls.config,cls.values)
        cls.rows,cls.ladder,cls.summary=cls.model.solve()

    def test_imposed_pressure_and_exact_main_pressure_event(self):
        before=[r for r in self.rows if r['time_s']<.02]
        np.testing.assert_allclose([r['torch_pressure_psia'] for r in before],150.,atol=1e-10)
        self.assertTrue(all(r['main_chamber_backpressure_psia']==14.8 for r in before))
        released=next(r for r in self.rows if r['time_s']==.02)
        self.assertEqual(released['pressure_imposed'],0)
        self.assertAlmostEqual(released['torch_pressure_psia'],150.)
        self.assertEqual(released['main_chamber_backpressure_psia'],300.)
        self.assertLess(released['throat_mdot_kg_s'],0)

    def test_backflow_and_mass_conservation(self):
        self.assertGreater(self.summary['main_chamber_backflow_kg'],0.)
        self.assertLess(self.summary['maximum_mass_balance_error_kg'],1e-9)
        active=[r for r in self.rows if r['time_s']>=.02]
        self.assertTrue(all(r['torch_gas_mass_kg']>0 and r['torch_temperature_k']>0 for r in active))
        self.assertGreater(self.summary['peak_purge_mdot_kg_s'],0.)
        steady=min(self.rows,key=lambda r:abs(r['time_s']-.1))
        self.assertGreater(steady['torch_pressure_psia'],300.)
        self.assertLess(steady['torch_pressure_psia'],301.)

    def test_steady_throat_balance_and_restriction_area_effect(self):
        equilibrium=self.model.fully_open_equilibrium(300.*psi)
        self.assertGreater(equilibrium['torch_over_main_pressure_psi'],0.)
        self.assertLess(equilibrium['torch_over_main_pressure_psi'],1.)
        self.assertFalse(equilibrium['throat_choked'])
        pressure=equilibrium['torch_pressure_psia']*psi
        supply=self.model.inlet(.1,pressure)[0]
        self.assertAlmostEqual(supply,self.model.throat(pressure,self.model.supply_temperature,300.*psi),places=10)
        smaller=TorchPurgeModel(self.config,self.values)
        smaller.throat_area *= .5
        self.assertGreater(smaller.fully_open_equilibrium(300.*psi)['torch_over_main_pressure_psi'],
                           equilibrium['torch_over_main_pressure_psi'])

    def test_timestep_refinement_preserves_purge_and_pressure_results(self):
        refined=copy.deepcopy(self.config)
        refined['torch']['standalone_purge']['maximum_step_s']=.0025
        refined['torch']['standalone_purge']['relative_tolerance']=1e-7
        rows,_,summary=TorchPurgeModel(refined,self.values).solve()
        self.assertAlmostEqual(summary['peak_released_torch_pressure_psia'],
                               self.summary['peak_released_torch_pressure_psia'],delta=.001)
        self.assertAlmostEqual(summary['total_purge_admitted_kg'],
                               self.summary['total_purge_admitted_kg'],delta=1e-8)
        self.assertAlmostEqual(summary['main_chamber_backflow_kg'],
                               self.summary['main_chamber_backflow_kg'],delta=1e-8)

    def test_ladder_contains_shared_fuel_injector_and_main_backpressure(self):
        self.assertEqual(self.ladder.shape[1],len(self.model.ladder_labels))
        self.assertIn('FIORF Outlet',self.model.ladder_labels)
        self.assertEqual(self.model.ladder_labels[-2:],[self.config['torch']['network']['chamber']['name'],'Main Chamber'])

    def test_missing_active_purge_area_is_rejected(self):
        config=copy.deepcopy(self.config)
        config['torch']['network']['purge_leg']['elements'][0]['cda_mm2']=None
        with self.assertRaisesRegex(ValueError,'CdA'):
            TorchPurgeModel(config,self.values)

    def test_purge_before_release_is_counted_without_changing_imposed_pressure(self):
        config=copy.deepcopy(self.config)
        valve=config['torch']['network']['purge_leg']['elements'][0]
        valve['command_events_s']=[[0.,1.],[.18,0.]]
        rows,_,summary=TorchPurgeModel(config,self.values).solve()
        self.assertGreater(summary['imposed_phase_purge_admitted_kg'],0.)
        self.assertGreater(summary['total_purge_admitted_kg'],summary['imposed_phase_purge_admitted_kg'])
        np.testing.assert_allclose([r['torch_pressure_psia'] for r in rows if r['time_s']<.02],150.)
        self.assertLess(summary['maximum_mass_balance_error_kg'],1e-9)


class CentralSequenceTests(unittest.TestCase):
    def setUp(self):
        self.config,self.values,_=load_config(ROOT/'hopper_config.json')

    def test_single_event_controls_valves_main_pressure_and_release_without_mutation(self):
        original=copy.deepcopy(self.config)
        model=TorchPurgeModel(self.config,self.values)
        self.assertEqual(self.config,original)
        self.assertEqual(model.release,2.)
        self.assertEqual(model.final_time,7.5)
        self.assertEqual(model.main_pressure(1.9)/psi,14.8)
        self.assertEqual(model.main_pressure(2.)/psi,300.)
        self.assertAlmostEqual(valve_position(model.valves['TFPV'],2.04),.5)
        self.assertAlmostEqual(valve_position(model.valves['FIMV'],2.04),.5)
        self.assertAlmostEqual(valve_position(model.valves['OIMV'],2.04),.6)
        event=self.config['torch']['standalone_purge']['sequence']['events'][2]
        event['time_s']=3.
        changed=TorchPurgeModel(self.config,self.values)
        self.assertEqual(changed.release,3.)
        self.assertEqual(changed.main_pressure(2.)/psi,14.8)
        self.assertEqual(valve_position(changed.valves['TFPV'],2.5),0.)

    def test_ambiguous_or_incomplete_sequences_fail(self):
        mutations=[
            lambda c:c['torch']['network']['purge_leg']['elements'][0].update(command_events_s=[[0.,0.]]),
            lambda c:c['torch']['standalone_purge'].update(final_time_s=7.5),
            lambda c:c['torch']['standalone_purge']['sequence']['events'][2]['valves'].update(TYPO=1),
            lambda c:c['torch']['standalone_purge']['sequence']['events'][2].update(time_s=.1),
            lambda c:c['torch']['standalone_purge']['sequence']['events'][2].pop('release_torch_pressure'),
            lambda c:c['torch']['standalone_purge']['sequence']['events'][0]['valves'].pop('TFPV'),
            lambda c:c['torch']['standalone_purge']['sequence']['events'][2].update(release_torch_pressure=False),
        ]
        for mutation in mutations:
            config=copy.deepcopy(self.config)
            mutation(config)
            with self.subTest(mutation=mutation),self.assertRaises(ValueError):
                TorchPurgeModel(config,self.values)


class TorchConnectionTests(unittest.TestCase):
    def test_purge_reaches_shared_injector_through_the_downstream_tee(self):
        config,_,_=load_config(ROOT/'hopper_config.json')
        rows=torch_connections(config['torch']['network'])
        edges={(r['from'],r['to']):r['solved_in_purge'] for r in rows}
        for edge in [('FIPORF','FIORF Outlet'),('FIORF Outlet','Torch Fuel Injector'),
                     ('Torch Fuel Injector','Torch Chamber'),('Torch Chamber','Torch Throat'),
                     ('Torch Throat','Main Chamber')]:
            self.assertTrue(edges[edge])
        self.assertFalse(edges[('FIORF','FIORF Outlet')])
        self.assertFalse(edges[('OTHX','Torch Ox Injector')])
        self.assertNotIn(('FIPORF','FIORF'),edges)
        self.assertNotIn(('FIPORF','Torch Chamber'),edges)

    def test_connection_export_follows_changed_tee_instead_of_hardcoding_fiorf(self):
        config,_,_=load_config(ROOT/'hopper_config.json')
        net=copy.deepcopy(config['torch']['network'])
        net['purge_leg']['sink']=net['fuel_leg']['elements'][1]['outlet_node']
        edges={(r['from'],r['to']) for r in torch_connections(net)}
        tee=net['purge_leg']['sink']
        self.assertIn(('FIPORF',tee),edges)
        self.assertIn((tee,'FIORF'),edges)


if __name__=='__main__':
    unittest.main()
