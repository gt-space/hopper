"""Independent cold-gas torch purge: quasi-steady restrictions + finite gas volume.

The first phase prescribes chamber pressure. Thereafter a well-mixed, adiabatic,
ideal-gas chamber integrates mass and internal energy. The downstream reservoir
can backflow through the geometric throat. Propellant flow/combustion, liquid
displacement and gas storage in the lines are deliberately outside this model.
"""
import math

import numpy as np
from scipy.constants import inch, liter, psi
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from thermoprop import Fluid

from Backend.Hopper_Torch_Config import validate_command_events, resolve_purge_sequence


def command_at(events, time):
    """Right-continuous command; before the first event the valve is closed."""
    index = np.searchsorted([row[0] for row in events], time, side='right') - 1
    return float(events[index][1]) if index >= 0 else 0.


def valve_position(spec, time):
    """Exact rate-limited response to step commands, including partial reversals."""
    events = spec['command_events_s']
    position, last_time, target = 0., 0., 0.

    def travel(position, target, elapsed):
        stroke = spec.get('open_time_s', 0.) if target >= position else spec.get('close_time_s', 0.)
        if stroke == 0:
            return target
        delta = max(0., elapsed) / stroke
        return min(target, position + delta) if target >= position else max(target, position - delta)

    for event_time, new_target in events:
        if event_time > time:
            break
        position = travel(position, target, event_time - last_time)
        target, last_time = new_target, event_time
        if spec.get('open_time_s' if target >= position else 'close_time_s', 0.) == 0:
            position = target
    return float(travel(position, target, time - last_time))


def orifice_flow(area, upstream, downstream, temperature, gamma, gas_constant):
    """Forward ideal-gas orifice flow (the FullFlow choked/subcritical equations)."""
    if area <= 0 or upstream <= downstream:
        return 0.
    ratio = max(0., downstream / upstream)
    critical = (2. / (gamma + 1.)) ** (gamma / (gamma - 1.))
    if ratio <= critical:
        factor = math.sqrt(gamma / (gas_constant * temperature)
                           * (2. / (gamma + 1.)) ** ((gamma + 1.) / (gamma - 1.)))
    else:
        # expm1 avoids cancellation near equal pressures.
        log_ratio = math.log(ratio)
        difference = math.exp(2./gamma*log_ratio) * -math.expm1((gamma-1.)/gamma*log_ratio)
        factor = math.sqrt(max(0., 2.*gamma/(gas_constant*temperature*(gamma-1.))*difference))
    return area * upstream * factor


def gas_chain(areas, upstream, downstream, temperature, gamma, gas_constant):
    """Solve common flow and each gas restriction's outlet pressure; no drop allocation.

    Nodes have no gas storage. A closed valve/check separates upstream and
    downstream equilibrated portions; multiply isolated trapped volumes are not
    represented. Backflow through the supply path is blocked by its check valve.
    """
    closed = next((i for i,a in enumerate(areas) if a <= 0), None)
    if closed is not None:
        return 0., [upstream]*(closed+1) + [downstream]*(len(areas)-closed)
    if downstream >= upstream:
        return 0., [upstream] + [downstream]*len(areas)
    critical = (2./(gamma+1.)) ** (gamma/(gamma-1.))
    sonic_factor = math.sqrt(gamma/(gas_constant*temperature)
                            * (2./(gamma+1.)) ** ((gamma+1.)/(gamma-1.)))

    def required_upstream(flow, outlet, area):
        sonic_pressure = flow / (area * sonic_factor)
        if sonic_pressure >= outlet/critical:
            return sonic_pressure
        if flow <= 0:
            return outlet
        return brentq(lambda p: orifice_flow(area,p,outlet,temperature,gamma,gas_constant)-flow,
                      outlet, outlet/critical, xtol=1e-5, rtol=1e-12)

    def pressures(flow):
        nodes = [downstream]
        for area in reversed(areas):
            nodes.append(required_upstream(flow,nodes[-1],area))
        return nodes[::-1]

    # Bracket slightly beyond the limiting choked capacity to tolerate roundoff
    # when the first restriction is choked and the endpoint residual is ~0.
    upper = min(areas)*upstream*sonic_factor*(1.+1e-8)
    if upper < 1e-14:
        cut = int(np.argmin(areas))
        return 0., [upstream]*(cut+1) + [downstream]*(len(areas)-cut)
    flow = brentq(lambda mass: pressures(mass)[0]-upstream,0.,upper,xtol=1e-13,rtol=1e-10)
    nodes = pressures(flow)
    nodes[0] = upstream
    return flow, nodes


def _positive(value, name):
    if isinstance(value,bool) or not isinstance(value,(float,int)) or not math.isfinite(value) or value <= 0:
        raise ValueError(f'{name} must be a positive finite number; fill the torch config first.')
    return float(value)


class TorchPurgeModel:
    def __init__(self, config, values):
        self.config, self.values = config, values
        self.torch = config['torch']
        self.network = self.torch.get('network')
        if self.network is None:
            raise ValueError('Independent torch purge needs torch.network, not legacy torch.branches.')
        self.settings, self.network = resolve_purge_sequence(self.torch.get('standalone_purge'), self.network)
        if not isinstance(self.settings,dict):
            raise ValueError('Add torch.standalone_purge settings; see README.')
        s, net = self.settings, self.network
        required = {'imposed_pressure_until_s','main_chamber_pressure_time_s_psia','main_chamber_gas_temperature_k',
                    'purge_supply_pressure_psia','purge_supply_temperature_k','gamma','final_time_s','maximum_step_s',
                    'save_interval_s','relative_tolerance','absolute_tolerance','near_zero_pressure_drop_pa',
                    'ladder_times_s','output_directory'}
        extra = s.keys()-required-{'note'}
        if required-s.keys() or extra:
            raise ValueError(f'standalone_purge: missing {sorted(required-s.keys())}, unknown {sorted(extra)}.')
        for key in required-{'purge_supply_pressure_psia','main_chamber_pressure_time_s_psia','ladder_times_s','output_directory','imposed_pressure_until_s'}:
            _positive(s[key], 'standalone_purge.'+key)
        if not isinstance(s['imposed_pressure_until_s'],(int,float)) or not math.isfinite(s['imposed_pressure_until_s']) or s['imposed_pressure_until_s'] < 0:
            raise ValueError('imposed_pressure_until_s must be finite and nonnegative.')
        if s['gamma'] <= 1 or s['imposed_pressure_until_s'] >= s['final_time_s']:
            raise ValueError('Need gamma > 1 and final_time_s after imposed_pressure_until_s.')
        self.release = s['imposed_pressure_until_s']
        self.final_time = s['final_time_s']
        if not isinstance(s['ladder_times_s'],list) or not s['ladder_times_s'] or any(
                not isinstance(t,(int,float)) or not math.isfinite(t) or not 0 <= t <= self.final_time for t in s['ladder_times_s']):
            raise ValueError('ladder_times_s must contain finite times within the run.')
        if not isinstance(s['output_directory'],str) or not s['output_directory'].strip():
            raise ValueError('standalone purge output_directory must be a nonempty path.')
        self.volume = _positive(net['chamber']['volume_l'],'torch chamber volume_l')*liter
        self.initial_pressure = _positive(net['chamber']['initial_pressure_psia'],'torch initial pressure')*psi
        self.initial_temperature = _positive(net['chamber']['initial_temperature_k'],'torch initial temperature')
        self.throat_area = _positive(net['throat']['area_in2'],'torch throat area_in2')*inch**2
        coefficient = net['throat'].get('discharge_coefficient')
        if coefficient is not None and coefficient != 1:
            raise ValueError('This purge model uses geometric throat area without an extra discharge coefficient.')
        self.main_events = s['main_chamber_pressure_time_s_psia']
        if not isinstance(self.main_events,list) or not self.main_events:
            raise ValueError('Main chamber pressure schedule must be a nonempty list.')
        previous = -1.
        for row in self.main_events:
            if (not isinstance(row,list) or len(row)!=2 or isinstance(row[0],bool)
                    or not isinstance(row[0],(int,float)) or not math.isfinite(row[0]) or row[0] <= previous):
                raise ValueError('Main chamber pressure times must be finite and strictly increasing.')
            _positive(row[1],'main chamber pressure')
            previous = row[0]
        if self.main_events[0][0] != 0:
            raise ValueError('Main chamber pressure schedule must start at zero.')
        if net['throat']['sink'] != 'Main Chamber':
            raise ValueError('Standalone purge currently requires the throat sink to be Main Chamber.')
        self.valves = {p['name']:p for key in ('fuel_leg','oxidizer_leg','purge_leg')
                       for p in net[key]['elements'] if p['type']=='Valve'}
        for spec in self.valves.values():
            if 'command_events_s' not in spec:
                raise ValueError(f'{spec["name"]}: add command_events_s for the standalone purge run.')
            validate_command_events(spec['command_events_s'])
        self.parts = list(net['purge_leg']['elements'])
        fuel = net['fuel_leg']['elements']
        tee = next(i for i,p in enumerate(fuel) if p['outlet_node']==net['purge_leg']['sink'])
        self.parts += fuel[tee+1:]
        if not any(p['type']=='CheckValve' for p in self.parts):
            raise ValueError('The purge supply path requires a check valve to block reverse flow.')
        if any(p['type'] not in ('Valve','CheckValve','Orifice','Injector') for p in self.parts):
            raise ValueError('Purge path supports valves, checks, orifices and injectors; no heated/multiphase tube.')
        self.areas = [_positive(p['cda_mm2'],p['name']+' CdA')*1e-6 for p in self.parts]
        self.supply_temperature = s['purge_supply_temperature_k']
        regulator = next((r for r in values['REGULATORS'] if r.outlet_node==net['purge_leg']['source']),None)
        if regulator is None:
            raise ValueError('Purge source must name a configured regulator outlet.')
        automatic_pressure = min(values['PRESSURANT_INITIAL_PRESSURE'], values['AMBIENT_PRESSURE']+regulator.set_pressure)
        self.supply_pressure = (automatic_pressure if s['purge_supply_pressure_psia'] is None
                                else _positive(s['purge_supply_pressure_psia'],'purge supply pressure')*psi)
        self.gamma = s['gamma']
        self.gas_constant = float(Fluid(values['PRESSURANT'],pressure=self.initial_pressure,
                                        temperature=self.initial_temperature).gas_constant)
        self.cv = self.gas_constant/(self.gamma-1.)
        self.cp = self.gamma*self.cv
        self.m0 = self.initial_pressure*self.volume/(self.gas_constant*self.initial_temperature)
        self.e0 = self.m0*self.cv*self.initial_temperature
        self.ladder_labels = [net['purge_leg']['source']] + [p['outlet_node'] for p in self.parts]
        self.ladder_labels[len(net['purge_leg']['elements'])] = net['purge_leg']['sink']
        self.ladder_labels[-1] = net['chamber']['name']
        self.ladder_labels.append(net['throat']['sink'])

    def main_pressure(self,time):
        events = self.main_events
        k = np.searchsorted([p[0] for p in events],time,side='right')-1
        return events[max(0,k)][1]*psi

    def sequence_rows(self):
        """One resolved command timeline for checking and exporting a run."""
        times = {0.,self.release,self.final_time}
        times.update(t for t,_ in self.main_events if t <= self.final_time)
        for spec in self.valves.values():
            times.update(t for t,_ in spec['command_events_s'] if t <= self.final_time)
        return [{'time_s':t,**{name:command_at(spec['command_events_s'],t) for name,spec in self.valves.items()},
                 'main_chamber_pressure_psia':self.main_pressure(t)/psi,
                 'torch_pressure_mode':'imposed initial pressure' if t < self.release else 'solve mass/energy',
                 'end_of_run':t == self.final_time} for t in sorted(times)]

    def inlet(self,time,pressure):
        areas = [a*valve_position(p,time) if p['type']=='Valve' else a for a,p in zip(self.areas,self.parts)]
        return gas_chain(areas,self.supply_pressure,pressure,self.supply_temperature,self.gamma,self.gas_constant)

    def throat(self,pressure,temperature,backpressure):
        upstream, downstream = max(pressure,backpressure), min(pressure,backpressure)
        source_temperature = temperature if pressure >= backpressure else self.settings['main_chamber_gas_temperature_k']
        drop = upstream-downstream
        threshold = self.settings['near_zero_pressure_drop_pa']
        if drop < threshold:
            flow = orifice_flow(self.throat_area,downstream+threshold,downstream,source_temperature,self.gamma,self.gas_constant)*drop/threshold
        else:
            flow = orifice_flow(self.throat_area,upstream,downstream,source_temperature,self.gamma,self.gas_constant)
        return flow if pressure >= backpressure else -flow

    def fully_open_equilibrium(self, backpressure):
        """Cold, adiabatic steady purge: inlet flow equals outward throat flow.

        With only one incoming gas stream, steady energy balance also gives
        chamber temperature = supply temperature. This is a diagnostic operating
        point, not a prescribed transient pressure or a design acceptance test.
        """
        def inlet(pressure):
            return gas_chain(self.areas,self.supply_pressure,pressure,self.supply_temperature,
                             self.gamma,self.gas_constant)[0]
        if self.supply_pressure <= backpressure:
            pressure, flow = backpressure, 0.
        else:
            pressure = brentq(lambda p: inlet(p)-self.throat(p,self.supply_temperature,backpressure),
                              backpressure,self.supply_pressure,xtol=1e-5)
            flow = inlet(pressure)
        return {'main_chamber_pressure_psia':backpressure/psi,'torch_pressure_psia':pressure/psi,
                'torch_over_main_pressure_psi':(pressure-backpressure)/psi,'purge_mdot_kg_s':flow,
                'throat_back_to_upstream_pressure_ratio':backpressure/pressure,
                'throat_choked':backpressure/pressure <= (2./(self.gamma+1.))**(self.gamma/(self.gamma-1.))}

    def rhs(self,time,state,backpressure):
        # Log states keep Newton/Jacobian trial mass and energy positive in this
        # very small (milliliter-scale) chamber, without clipping accepted states.
        mass, energy = math.exp(state[0])*self.m0, math.exp(state[1])*self.e0
        temperature = energy/(mass*self.cv)
        pressure = (self.gamma-1.)*energy/self.volume
        inlet,_ = self.inlet(time,pressure)
        throat = self.throat(pressure,temperature,backpressure)
        incoming = max(-throat,0.)
        outgoing = max(throat,0.)
        heat = (inlet*self.cp*self.supply_temperature
                + incoming*self.cp*self.settings['main_chamber_gas_temperature_k'] - outgoing*self.cp*temperature)
        # Last three states integrate admitted purge, outward throat flow and backflow.
        return [(inlet-throat)/mass,heat/energy,inlet/self.m0,outgoing/self.m0,incoming/self.m0]

    def breakpoints(self, include_imposed=False):
        points = {0.,self.release,self.final_time}
        points.update(t for t,_ in self.main_events)
        for spec in self.valves.values():
            points.update(t for t,_ in spec['command_events_s'])
            for t,_ in spec['command_events_s']:
                points.update((t+spec.get('open_time_s',0.),t+spec.get('close_time_s',0.)))
        return sorted(t for t in points if (0. if include_imposed else self.release) <= t <= self.final_time)

    def solve(self):
        s = self.settings
        time = np.arange(0.,self.final_time+s['save_interval_s']*.1,s['save_interval_s'])
        # Resolve the rapid finite-volume response to pressure steps in the exported data.
        extra = [event+offset for event,_ in self.main_events if event >= self.release
                 for offset in (0.,1e-6,1e-5,1e-4,5e-4,1e-3)]
        time = np.unique(np.r_[time,[0.,self.release,self.final_time],extra,self.breakpoints(),s['ladder_times_s']])
        time = time[(time >= 0)&(time <= self.final_time)]
        state = np.zeros(5)
        history = np.tile(state,(len(time),1))
        # While pressure is prescribed, account for any scheduled early purge
        # without pretending to solve the unmodelled combustion/outlet balance.
        imposed_breaks = [t for t in self.breakpoints(include_imposed=True) if t <= self.release]
        for start,end in zip(imposed_breaks,imposed_breaks[1:]):
            solution = solve_ivp(lambda t,y:[self.inlet(t,self.initial_pressure)[0]/self.m0],
                                 (start,end),[state[2]],rtol=s['relative_tolerance'],
                                 atol=s['absolute_tolerance'],max_step=s['maximum_step_s'],dense_output=True)
            if not solution.success:
                raise RuntimeError('Failed to integrate purge admission during the imposed-pressure phase.')
            mask = (time >= start)&(time <= end)
            history[mask,2] = solution.sol(time[mask])[0]
            state[2] = solution.y[0,-1]
        imposed_purge = state[2]
        for start,end in zip(self.breakpoints(),self.breakpoints()[1:]):
            backpressure = self.main_pressure(start)
            solution = solve_ivp(lambda t,y:self.rhs(t,y,backpressure),(start,end),state,method='BDF',
                                 rtol=s['relative_tolerance'],atol=s['absolute_tolerance'],
                                 max_step=s['maximum_step_s'],dense_output=True)
            if not solution.success:
                raise RuntimeError(f'Torch purge solve failed at {solution.t[-1]:g} s: {solution.message}')
            mask = (time >= start)&(time <= end)
            history[mask] = solution.sol(time[mask]).T
            state = solution.y[:,-1]
        rows, pressures = [],[]
        for t,y in zip(time,history):
            mass,energy = math.exp(y[0])*self.m0,math.exp(y[1])*self.e0
            pressure=(self.gamma-1.)*energy/self.volume
            temperature=energy/(mass*self.cv)
            inlet,nodes = self.inlet(t,pressure)
            imposed=t < self.release
            backpressure=self.main_pressure(t)
            throat=self.throat(pressure,temperature,backpressure)
            row = {'time_s':float(t),'pressure_imposed':int(imposed),'torch_pressure_psia':pressure/psi,
                   'main_chamber_backpressure_psia':backpressure/psi,'purge_supply_pressure_psia':self.supply_pressure/psi,
                   'torch_over_main_pressure_psi':(pressure-backpressure)/psi,
                   'torch_temperature_k':temperature,'purge_mdot_kg_s':inlet,
                   'throat_mdot_kg_s':float('nan') if imposed else throat,
                   'torch_gas_mass_kg':float('nan') if imposed else mass,
                   'purge_admitted_kg':y[2]*self.m0,'throat_outward_kg':y[3]*self.m0,
                   'main_chamber_backflow_kg':y[4]*self.m0}
            for name,spec in self.valves.items():
                row[name+'_command']=command_at(spec['command_events_s'],t)
                row[name+'_opening']=valve_position(spec,t)
            rows.append(row)
            pressures.append(nodes+[backpressure])
        active=history[time>=self.release]
        residual=(np.exp(active[:,0])-1.-(active[:,2]-imposed_purge)+active[:,3]-active[:,4])*self.m0
        summary = {'model':'Cold ideal-gas purge; no combustion, propellant flow or liquid displacement',
                   'imposed_pressure_until_s':self.release,'initial_torch_pressure_psia':self.initial_pressure/psi,
                   'purge_source':self.network['purge_leg']['source'],'purge_supply_pressure_psia':self.supply_pressure/psi,
                   'initial_temperature_k':self.initial_temperature,'main_reservoir_temperature_k':s['main_chamber_gas_temperature_k'],
                   'imposed_phase_purge_admitted_kg':float(imposed_purge*self.m0),
                   'total_purge_admitted_kg':float(state[2]*self.m0),'main_chamber_backflow_kg':float(state[4]*self.m0),
                   'maximum_mass_balance_error_kg':float(np.max(np.abs(residual))),
                   'peak_purge_mdot_kg_s':max(r['purge_mdot_kg_s'] for r in rows),
                   'minimum_throat_mdot_kg_s':min(r['throat_mdot_kg_s'] for r in rows if not r['pressure_imposed']),
                   'peak_released_torch_pressure_psia':max(r['torch_pressure_psia'] for r in rows if not r['pressure_imposed']),
                   'fully_open_purge_equilibria':[self.fully_open_equilibrium(p*psi)
                                                  for p in sorted({self.main_pressure(self.release)/psi}
                                                  | {p for t,p in self.main_events if self.release <= t <= self.final_time})],
                   'assumptions':['Constant regulated pressure source; bottle depletion/regulator droop excluded.',
                                  'Both reservoirs and the torch use the same ideal nitrogen model with constant gamma.',
                                  'Mass and energy integrate only after the imposed-pressure phase.',
                                  'Fuel and oxidizer flow, closing-valve residual propellant and combustion excluded.',
                                  'No trapped gas/remaining liquid inventory in the zero-volume purge/fuel line.',
                                  'Throat is bidirectional; forward sign is torch to main chamber.',
                                  'Near-zero throat pressure drop is linearly regularized below the configured threshold.']}
        return rows,np.asarray(pressures)/psi,summary
