"""Validated torch topology inputs, separate from the legacy lumped branches.

This describes the requested network; it does not assemble FullFlow components.
Unset hardware measurements are represented by None, never invented geometry.
"""
import inspect
import math
import copy

import numpy as np

from Backend import Hopper_Elements as elements
from Backend.Hopper_Config_Fields import UNITS


def _keys(value, required, optional=(), where='torch.network'):
    if not isinstance(value, dict):
        raise ValueError(f'{where}: expected an object.')
    missing = set(required) - value.keys()
    extra = value.keys() - set(required) - set(optional)
    if missing or extra:
        raise ValueError(f'{where}: missing {sorted(missing)}, unknown {sorted(extra)}.')


def _positive_or_unset(value, where):
    if value is not None and (isinstance(value, bool) or not isinstance(value, (int, float))
                              or not math.isfinite(value) or value <= 0):
        raise ValueError(f'{where}: must be positive or null (not yet specified).')


def validate_time_curve(points):
    if not isinstance(points, list) or not points:
        raise ValueError('OTHX time_s_delta_t_k requires at least one [time_s, delta_t_k] point.')
    previous = -math.inf
    for row in points:
        if not isinstance(row, list) or len(row) != 2:
            raise ValueError('OTHX time_s_delta_t_k points must be [time_s, delta_t_k].')
        if any(isinstance(x, bool) or not isinstance(x, (int, float)) or not math.isfinite(x) for x in row):
            raise ValueError('OTHX time/dT values must be finite numbers.')
        t, _ = row
        if t < 0 or t <= previous:
            raise ValueError('OTHX times must be nonnegative and strictly increasing.')
        previous = t


def validate_command_events(events):
    """Step commands [simulation seconds, opening 0..1], followed by physical travel."""
    validate_time_curve(events)
    if any(not 0 <= opening <= 1 for _, opening in events):
        raise ValueError('Valve command_events_s openings must be between zero and one.')


def resolve_purge_sequence(settings, network):
    """Expand one human-readable timeline; retain old per-valve archives for replay.

    The returned objects are private copies. Original JSON remains the single
    source for archiving; expanded valve schedules cannot shadow user inputs.
    """
    settings, network = copy.deepcopy(settings), copy.deepcopy(network)
    if not isinstance(settings,dict) or 'sequence' not in settings:
        return settings, network
    legacy = {'imposed_pressure_until_s','main_chamber_pressure_time_s_psia','final_time_s'}
    if legacy & settings.keys():
        raise ValueError('Use standalone_purge.sequence only; remove legacy pressure/timing fields.')
    sequence = settings.pop('sequence')
    _keys(sequence, ('end_time_s','events'), ('note',), 'standalone_purge.sequence')
    end = sequence['end_time_s']
    _positive_or_unset(end, 'sequence.end_time_s')
    if end is None:
        raise ValueError('sequence.end_time_s must be specified.')
    events = sequence['events']
    if not isinstance(events,list) or not events:
        raise ValueError('sequence.events must be a nonempty list.')
    valves = {p['name']:p for key in ('fuel_leg','oxidizer_leg','purge_leg')
              for p in network[key]['elements'] if p['type']=='Valve'}
    if any('command_events_s' in p or p.get('schedule') is not None for p in valves.values()):
        raise ValueError('Central sequence owns all torch valve commands; remove component command_events_s/schedule.')
    commands = {name:[] for name in valves}
    main, release, previous = [], [], -1.
    for event in events:
        _keys(event, ('time_s',), ('valves','main_chamber_pressure_psia','release_torch_pressure','note'), 'sequence event')
        t = event['time_s']
        if (isinstance(t,bool) or not isinstance(t,(int,float)) or not math.isfinite(t)
                or t < 0 or t <= previous or t > end):
            raise ValueError('sequence event times must increase, be nonnegative, and not exceed end_time_s.')
        previous = t
        changes = event.get('valves',{})
        if not isinstance(changes,dict) or changes.keys()-valves.keys():
            raise ValueError('sequence valves must name configured torch valves only.')
        for name,opening in changes.items():
            validate_command_events([[t,opening]])
            commands[name].append([t,opening])
        if 'main_chamber_pressure_psia' in event:
            pressure = event['main_chamber_pressure_psia']
            _positive_or_unset(pressure,'sequence main_chamber_pressure_psia')
            if pressure is None:
                raise ValueError('sequence main pressure must be specified.')
            main.append([t,pressure])
        if 'release_torch_pressure' in event:
            if event['release_torch_pressure'] is not True:
                raise ValueError('release_torch_pressure must be true at exactly one event.')
            release.append(t)
    if any(not rows or rows[0][0] != 0 for rows in commands.values()):
        raise ValueError('Initialize every torch valve at time_s=0 in sequence.events.')
    if not main or main[0][0] != 0:
        raise ValueError('Initialize main_chamber_pressure_psia at time_s=0 in sequence.events.')
    if len(release) != 1 or release[0] >= end:
        raise ValueError('sequence needs exactly one release_torch_pressure=true event before end_time_s.')
    for name,part in valves.items():
        part['command_events_s'] = commands[name]
    settings.update(imposed_pressure_until_s=release[0],main_chamber_pressure_time_s_psia=main,final_time_s=end)
    return settings, network


def heat_exchanger_outlet_temperature(spec, inlet_temperature, simulation_time, torch_open_time):
    """Configured Tout=Tin+dT(t), in K/s; linear interpolation, held endpoints.

    Torch-relative curves apply no heating before the torch opening command.
    This helper is for the forthcoming torch builder, not the current main solver.
    """
    heating = spec['temperature_rise']
    validate_time_curve(heating['time_s_delta_t_k'])
    if heating['time_basis'] == 'torch_open':
        t = simulation_time - torch_open_time
        if t < 0:
            return float(inlet_temperature)
    elif heating['time_basis'] == 'simulation':
        t = simulation_time
    else:
        raise ValueError('OTHX time_basis must be torch_open or simulation.')
    points = heating['time_s_delta_t_k']
    rise = np.interp(t, [row[0] for row in points], [row[1] for row in points])
    return float(inlet_temperature + rise)


def _component(spec):
    if not isinstance(spec, dict):
        raise ValueError('torch component must be an object.')
    kind = spec.get('type')
    if kind not in ('Valve', 'CheckValve', 'Orifice', 'Injector', 'HeatedLine'):
        raise ValueError(f'torch.network: unsupported component type {kind!r}.')
    cls = elements.Line if kind == 'HeatedLine' else getattr(elements, kind)
    parameters = inspect.signature(cls.__init__).parameters
    allowed = {'type', 'note'} | {UNITS.get(name, (name, 1))[0] for name in parameters if name != 'self'}
    if kind == 'HeatedLine':
        allowed.add('temperature_rise')
    if kind == 'Valve':
        allowed.add('command_events_s')
    _keys(spec, ('type', 'name'), allowed - {'type', 'name'}, 'torch component')
    if not isinstance(spec['name'], str) or not spec['name'].strip():
        raise ValueError('torch component needs a nonempty name.')
    spec.setdefault('outlet_node', spec['name'] + ' Outlet')
    if not isinstance(spec['outlet_node'], str) or not spec['outlet_node'].strip():
        raise ValueError(f'{spec["name"]}: outlet_node must be a nonempty name.')
    for key in ('cda_mm2', 'length_in', 'inner_diameter_in', 'outer_diameter_in', 'wall_thickness_in'):
        if key in spec:
            _positive_or_unset(spec[key], f'{spec["name"]}.{key}')
    for key in ('open_time_s', 'close_time_s', 'roughness_in'):
        if key in spec:
            value = spec[key]
            if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value < 0:
                raise ValueError(f'{spec["name"]}.{key}: must be finite and nonnegative.')
    if kind in ('Valve', 'CheckValve'):
        if 'command_events_s' in spec:
            validate_command_events(spec['command_events_s'])
            if spec.get('schedule') is not None:
                raise ValueError(f'{spec["name"]}: choose command_events_s or schedule, not both.')
        schedule = spec.get('schedule')
        if isinstance(schedule, str):
            if schedule not in ('torch', 'closed'):
                raise ValueError(f'{spec["name"]}: unknown valve schedule {schedule!r}.')
        elif schedule is not None:
            if (not isinstance(schedule, list) or len(schedule) != 2
                    or not schedule[0] or len(schedule[0]) != len(schedule[1])
                    or any(not math.isfinite(x) or x < 0 for x in schedule[0])
                    or any(b <= a for a, b in zip(schedule[0], schedule[0][1:]))
                    or any(not math.isfinite(x) or not 0 <= x <= 1 for x in schedule[1])):
                raise ValueError(f'{spec["name"]}: invalid valve schedule.')
    if kind == 'HeatedLine':
        _keys(spec, ('type', 'name', 'length_in', 'temperature_rise'), allowed - {'type', 'name', 'length_in', 'temperature_rise'}, 'OTHX')
        od, wall, diameter = (spec.get(k) for k in ('outer_diameter_in', 'wall_thickness_in', 'inner_diameter_in'))
        if od is not None and wall is not None:
            calculated = od - 2*wall
            if calculated <= 0 or (diameter is not None and not math.isclose(calculated, diameter, rel_tol=1e-6)):
                raise ValueError('OTHX: ID must agree with positive OD minus twice the wall.')
        heating = spec['temperature_rise']
        _keys(heating, ('time_basis', 'time_s_delta_t_k'), ('note',), 'OTHX.temperature_rise')
        if heating['time_basis'] not in ('torch_open', 'simulation'):
            raise ValueError('OTHX time_basis must be torch_open or simulation.')
        validate_time_curve(heating['time_s_delta_t_k'])


def validate_torch_network(network, sources):
    """Validate declared nodes/tee connectivity, retaining unmeasured null inputs."""
    _keys(network, ('fuel_leg', 'purge_leg', 'oxidizer_leg', 'chamber', 'throat'), ('note',))
    chamber, throat = network['chamber'], network['throat']
    _keys(chamber, ('name', 'volume_l', 'initial_pressure_psia', 'initial_temperature_k'), ('note',), 'torch chamber')
    _keys(throat, ('name', 'source', 'sink', 'area_in2'), ('note','discharge_coefficient'), 'torch throat')
    if not isinstance(chamber['name'], str) or not chamber['name'] or not isinstance(throat['name'], str) or not throat['name']:
        raise ValueError('Torch chamber/throat require nonempty names.')
    for key in ('volume_l', 'initial_pressure_psia', 'initial_temperature_k'):
        _positive_or_unset(chamber[key], 'torch.chamber.' + key)
    for key in ('area_in2', 'discharge_coefficient'):
        _positive_or_unset(throat.get(key), 'torch.throat.' + key)
    if throat['source'] != chamber['name'] or throat['sink'] not in sources:
        raise ValueError('Torch throat must connect the torch chamber to a known external sink.')
    names = [chamber['name'], throat['name']]
    outlets = set()
    for key in ('fuel_leg', 'purge_leg', 'oxidizer_leg'):
        leg = network[key]
        _keys(leg, ('source', 'sink', 'elements'), ('note',), key)
        if leg['source'] not in sources:
            raise ValueError(f'{key}: unknown source {leg["source"]!r}.')
        if not isinstance(leg['elements'], list) or not leg['elements']:
            raise ValueError(f'{key}: requires an ordered, nonempty elements list.')
        for part in leg['elements']:
            _component(part)
            names.append(part['name'])
            if part['outlet_node'] in outlets or part['outlet_node'] in sources or part['outlet_node'] == chamber['name']:
                raise ValueError(f'Duplicate/reserved torch outlet node {part["outlet_node"]!r}.')
            outlets.add(part['outlet_node'])
    if len(set(names)) != len(names):
        raise ValueError('Torch chamber, throat, and component names must be unique.')
    for key in ('fuel_leg', 'oxidizer_leg'):
        if network[key]['sink'] != chamber['name']:
            raise ValueError(f'{key} must terminate at the shared torch chamber.')
    fuel_outlets = {e['outlet_node'] for e in network['fuel_leg']['elements'][:-1]}
    if network['purge_leg']['sink'] not in fuel_outlets:
        raise ValueError('purge_leg must tee into an intermediate fuel-leg outlet, upstream of its injector.')
    return network
