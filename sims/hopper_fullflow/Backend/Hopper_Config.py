"""JSON configuration boundary. Physics continues to consume SI-valued inputs."""

import inspect
import json
import math
import os
from pathlib import Path

from scipy.constants import psi

from Backend import Hopper_Elements as elements
from Backend.Hopper_Inputs import LBF, load_thrust_curve, ramp_profile, table_value, validate_profile

# Unit conversions and solver bindings are code, not vehicle settings.
from Backend.Hopper_Config_Fields import GROUP_EXTRAS, GROUPS, UNITS
from Backend.Hopper_Torch_Config import validate_torch_network

ROOT = Path(__file__).resolve().parents[1]   # hopper_fullflow
SCHEMA_VERSION = 2
TYPES = {name: getattr(elements, name) for name in (
    'Tank', 'COPVSource', 'Line', 'GasLine', 'Orifice', 'Injector', 'Regen',
    'Valve', 'BangBangValve', 'CheckValve', 'ThrottleValve', 'CavitatingVenturi',
    'ReliefValve', 'Regulator')}
# Top-level sections besides GROUPS.
SECTIONS = {'schema_version', 'thrust_profile', 'bottles', 'tanks', 'flow_paths', 'regulators',
            'side_paths', 'torch', 'purge', 'not_modelled', 'monte_carlo'}
OPTIONAL_SECTIONS = {'description', 'notes'}


def _keys(value, allowed, where, required=()):
    if not isinstance(value, dict):
        raise ValueError(f'{where}: expected an object.')
    extra, missing = value.keys() - set(allowed), set(required) - value.keys()
    if extra or missing:
        raise ValueError(f'{where}: unknown fields {sorted(extra)}; missing fields {sorted(missing)}.')


def _unique_object(pairs):
    value = {}
    for key, item in pairs:
        if key in value:
            raise ValueError(f'Duplicate JSON field: {key}')
        value[key] = item
    return value


def _scaled(raw, scale):
    if scale == 1 or raw is None:
        return raw
    if isinstance(raw, list):
        return [item * scale for item in raw]
    return raw * scale


def build_component(spec, schedules=None):
    """A hardware object from its JSON spec. "note" (free text, e.g. "EDIT GUESS ...")
    is kept on the object for the component table and the parameter workbook."""
    kind = spec.get('type')
    if kind not in TYPES:
        raise ValueError(f'Unknown component type {kind!r}; choose from {list(TYPES)}.')
    cls = TYPES[kind]
    schedules = {'closed': elements.CLOSED, **(schedules or {})}
    parameters = inspect.signature(cls.__init__).parameters
    allowed = {'type', 'note'}
    kwargs = {}
    for name, parameter in parameters.items():
        if name == 'self':
            continue
        key, scale = UNITS.get(name, (name, 1))
        allowed.add(key)
        if key not in spec:
            if parameter.default is inspect.Parameter.empty:
                raise ValueError(f'{kind}: missing {key}.')
            spec[key] = parameter.default
        value = spec[key]
        if name == 'delta_t_curve':
            value = [(p * psi, t) for p, t in value]
        elif name == 'schedule' and isinstance(value, str):
            if value not in schedules:
                raise ValueError(f'{spec.get("name", kind)}: unknown schedule {value!r}; choose from {list(schedules)}.')
            value = schedules[value]
        elif value is not None and scale != 1:
            value *= scale
        kwargs[name] = value
    _keys(spec, allowed, spec.get('name', kind))
    for field in ('open_time', 'close_time', 'roughness'):
        if field in kwargs and kwargs[field] < 0:
            raise ValueError(f'{kwargs["name"]}: {field} cannot be negative.')
    for field in ('volume', 'pressure', 'liquid_mass', 'liquid_temperature', 'ullage_temperature',
                  'set_pressure', 'set_gauge_pressure', 'crack_pressure', 'cda', 'max_cda', 'throat_area'):
        if kwargs.get(field) is not None and kwargs[field] <= 0:
            raise ValueError(f'{kwargs["name"]}: {field} must be positive.')
    if 'collapse_factor' in kwargs and kwargs['collapse_factor'] < 1:
        raise ValueError(f'{kwargs["name"]}: collapse_factor below 1 would create gas from nowhere.')
    for field in ('upper_limit', 'lower_limit'):
        if field in kwargs and kwargs[field] < 0:
            raise ValueError(f'{kwargs["name"]}: deadband cannot be negative.')
    schedule = kwargs.get('schedule')
    if schedule is not None:
        times, openings = schedule
        if not times or len(times) != len(openings) or any(b <= a for a,b in zip(times,times[1:])) or any(not 0 <= x <= 1 for x in openings):
            raise ValueError(f'{kwargs["name"]}: invalid valve schedule [times_s, openings].')
    component = cls(**kwargs)
    component.note = spec.get('note', '')
    return component


def build_branch(spec, schedules=None):
    """A side path (vent, relief, purge, drain, torch feed) from its JSON spec."""
    _keys(spec, ('type', 'name', 'source', 'sink', 'into_sink', 'elements', 'note'),
          spec.get('name', 'Branch'), ('name', 'source', 'elements'))
    spec.setdefault('sink', 'ambient')
    spec.setdefault('into_sink', False)
    parts = [build_component(element, schedules) for element in spec['elements']]
    if not parts or any(not isinstance(part, (elements.GasLine, elements.Orifice)) for part in parts):
        raise ValueError(f'{spec["name"]}: a branch needs one or more GasLine / Orifice / Valve / relief elements.')
    branch = elements.Branch(spec['name'], spec['source'], parts, spec['sink'], spec['into_sink'])
    branch.note = spec.get('note', '')
    return branch


def efficiency(spec, base):
    mode = spec['mode']
    if mode == 'linear':
        _keys(spec, ('mode','low','high','pc_low_psia','pc_span_psi','csv'), 'efficiency', ('low','high','pc_low_psia','pc_span_psi'))
        if spec['pc_span_psi'] <= 0 or not all(0 < spec[k] <= 1 for k in ('low','high')):
            raise ValueError('Efficiency endpoints must be in (0,1], with positive pressure span.')
        def function(pc, mr):
            f = min(1., max(0., (pc / psi - spec['pc_low_psia']) / spec['pc_span_psi']))
            return spec['low'] + (spec['high'] - spec['low']) * f
    elif mode == 'polynomial':
        _keys(spec, ('mode','coefficients'), 'efficiency', ('coefficients',))
        # coefficients[i][j] multiplies (Pc in psia)^i * MR^j.
        coefficients = spec['coefficients']
        if not coefficients or any(not row for row in coefficients):
            raise ValueError('Efficiency polynomial requires a nonempty coefficient matrix.')
        def function(pc, mr):
            return sum(c * (pc / psi)**i * mr**j for i,row in enumerate(coefficients) for j,c in enumerate(row))
    elif mode == 'csv':
        _keys(spec, ('mode','csv'), 'efficiency', ('csv',))
        function = None
    else:
        raise ValueError('Efficiency mode must be linear, polynomial, or csv.')
    filename = str((base / spec['csv']).resolve()) if spec.get('csv') else None
    if mode == 'csv' and not filename:
        raise ValueError('CSV efficiency requires a filename.')
    if filename:
        table_value(filename, 100 * psi, 1.)  # fail early on invalid/missing tables
    return function, filename


def _thrust_profiles(profile, base):
    _keys(profile, ('selected','interpolation','csv_min_spacing_s','ramp','profiles'), 'thrust_profile',
          ('selected','interpolation','csv_min_spacing_s','ramp','profiles'))
    if profile['interpolation'] not in ('linear','step'):
        raise ValueError('Thrust interpolation must be linear or step.')
    ramp = profile['ramp']
    _keys(ramp, ('start_lbf','end_lbf','duration_s','points'), 'ramp', ('start_lbf','end_lbf','duration_s','points'))
    curves = {}
    for name, entry in profile['profiles'].items():
        if name == 'ramp':
            raise ValueError('The profile name ramp is reserved for the configured ramp.')
        if 'csv' in entry:
            _keys(entry, ('csv','note'), name, ('csv',))
            curves[name] = validate_profile(load_thrust_curve(base / entry['csv'], profile['csv_min_spacing_s']))
            continue
        _keys(entry, ('time_s_thrust_n','note'), name, ('time_s_thrust_n',))
        pairs = entry['time_s_thrust_n']
        if len(dict(pairs)) != len(pairs):
            raise ValueError(f'{name}: duplicate thrust times.')
        curves[name] = validate_profile(dict(pairs))
    curves['ramp'] = ramp_profile(ramp['start_lbf']*LBF, ramp['end_lbf']*LBF, ramp['duration_s'], ramp['points'])
    if profile['selected'] not in curves:
        raise ValueError('Selected thrust profile is not defined.')
    return curves


def _monte_carlo(mc, base):
    _keys(mc, ('runs','seed','workers','save_interval_s','keep_files','directory','thrust_control','dispersions','input_limits'),
          'monte_carlo', ('runs','seed','workers','save_interval_s','keep_files','directory','thrust_control','dispersions','input_limits'))
    if mc['thrust_control'] not in ('thrust','flow'):
        raise ValueError('monte_carlo.thrust_control must be thrust or flow.')
    dispersions = {}
    for group, entry in mc['dispersions'].items():
        _keys(entry, ('distribution','size','mode','inputs','how','note'), group, ('distribution','size','mode','inputs','how'))
        if entry['distribution'] not in ('uniform','normal') or entry['mode'] not in ('rel','abs') or entry['how'] not in ('each','shared'):
            raise ValueError(f'{group}: distribution must be uniform/normal, mode rel/abs, how each/shared.')
        dispersions[group] = (entry['distribution'], entry['size'], entry['mode'], entry['inputs'], entry['how'])
    return {'MC_RUNS': mc['runs'], 'MC_SEED': mc['seed'], 'MC_WORKERS': mc['workers'], 'MC_SAVE_DT': mc['save_interval_s'],
            'MC_KEEP_FILES': mc['keep_files'], 'MC_DIRECTORY': str((base / mc['directory']).resolve()),
            'MC_THRUST_CONTROL': mc['thrust_control'], 'MC_ALL_CDAS': 'ALL_CDAS', 'MC_DISPERSIONS': dispersions,
            'MC_INPUT_LIMITS': {name: tuple(limits) for name, limits in mc['input_limits'].items()}}


def _purge(spec, branches):
    """The post-shutdown purge sequence."""
    keys = ('sequence','valve_stroke_s','simulate_in_transient')
    _keys(spec, keys, 'purge', keys)
    by_name = {branch.name: branch for branch in branches}
    sequence = []
    if spec['valve_stroke_s'] <= 0:
        raise ValueError('purge.valve_stroke_s must be positive.')
    for leg in spec['sequence']:
        _keys(leg, ('branch','open_after_shutdown_s','duration_s','note'), 'purge.sequence',
              ('branch','open_after_shutdown_s','duration_s'))
        branch = by_name.get(leg['branch'])
        if branch is None:
            raise ValueError(f"purge.sequence: no side path called {leg['branch']!r}.")
        if not any(isinstance(e, elements.Valve) and not isinstance(e, elements.CheckValve) for e in branch.elements):
            raise ValueError(f"purge.sequence: {leg['branch']} has no valve to open.")
        if leg['duration_s'] <= spec['valve_stroke_s'] or leg['open_after_shutdown_s'] < 0:
            raise ValueError(f"purge.sequence: {leg['branch']} must open at or after shutdown and stay open "
                             "longer than the valve stroke.")
        sequence.append((leg['branch'], float(leg['open_after_shutdown_s']), float(leg['duration_s'])))
    if len({name for name, _, _ in sequence}) != len(sequence):
        raise ValueError('purge.sequence: each branch can appear once.')
    return {'PURGE_SEQUENCE': sequence, 'PURGE_VALVE_STROKE': spec['valve_stroke_s'],
            'PURGE_IN_TRANSIENT': spec['simulate_in_transient']}


TORCH_SEQUENCE_KEYS = ('oxidizer_open_s', 'fuel_open_s', 'close_after_ignition_s', 'purge_open_after_ignition_s',
                       'purge_close_after_main_purges_s')


def _torch_legs(network):
    """Element objects for the enabled torch network: {"fuel" | "oxidizer" | "purge":
    {"source", "elements", "stroke"}}. Valve schedules are set by the transient from
    torch.sequence; each valve's open_time_s / close_time_s in the network is its
    stroke (the schedule's opening and closing ramp). OTHX (HeatedLine) is built as a
    tube with its friction-equivalent CdA (no heating). The purge leg tees in
    upstream of the torch fuel injector; legs["purge"]["shared"] lists the fuel-leg
    elements after the tee, which the purge gas also flows through."""
    legs = {}
    for label, key in (('fuel', 'fuel_leg'), ('oxidizer', 'oxidizer_leg'), ('purge', 'purge_leg')):
        leg = network[key]
        parts, strokes = [], {}
        for raw in leg['elements']:
            spec = {k: v for k, v in raw.items() if k not in ('command_events_s', 'temperature_rise', 'outlet_node')}
            if spec['type'] == 'HeatedLine':
                spec['type'] = 'GasLine'
                spec.pop('height_change_in', None)
            if spec['type'] == 'Valve':
                strokes[spec['name']] = (spec.pop('open_time_s', 0.0) or 0.0, spec.pop('close_time_s', 0.0) or 0.0)
                spec['schedule'] = 'closed'
            if any(spec.get(k, 0) is None for k in ('cda_mm2', 'length_in')):
                raise ValueError(f'torch.enabled: {spec["name"]} needs every CdA / length set (no nulls).')
            parts.append(build_component(spec))
        legs[label] = {'source': leg['source'], 'elements': parts, 'stroke': strokes}
    fuel = network['fuel_leg']['elements']
    tee = next(i for i, e in enumerate(fuel) if e['outlet_node'] == network['purge_leg']['sink'])
    legs['purge']['shared'] = legs['fuel']['elements'][tee + 1:]
    if network['chamber']['volume_l'] is None or network['throat']['area_in2'] is None:
        raise ValueError('torch.enabled: the torch chamber volume and throat area must be set.')
    legs['chamber_name'] = network['chamber']['name']
    legs['throat_area'] = network['throat']['area_in2'] * 0.0254**2
    legs['throat_discharge_coefficient'] = network['throat'].get('discharge_coefficient') or 1.0
    return legs


def _torch_sequence(spec, values):
    """Flight torch timing (startup and the torch purge), simulation seconds."""
    if spec is None:
        raise ValueError(f'torch.enabled with a network needs torch.sequence {list(TORCH_SEQUENCE_KEYS)}.')
    _keys(spec, TORCH_SEQUENCE_KEYS + ('note',), 'torch.sequence', TORCH_SEQUENCE_KEYS)
    if any(isinstance(spec[k], bool) or not isinstance(spec[k], (int, float)) or spec[k] < 0 for k in TORCH_SEQUENCE_KEYS):
        raise ValueError('torch.sequence: times must be nonnegative numbers.')
    if max(spec['oxidizer_open_s'], spec['fuel_open_s']) >= values['PROFILE_START_TIME'] + spec['close_after_ignition_s']:
        raise ValueError('torch.sequence: the torch valves must open before they close.')
    return {k: float(spec[k]) for k in TORCH_SEQUENCE_KEYS}


def load_config(filename=None):
    path = Path(filename or os.environ.get('HOPPER_CONFIG', ROOT / 'hopper_config.json')).resolve()
    config = json.loads(path.read_text(encoding='utf-8-sig'), object_pairs_hook=_unique_object,
                        parse_constant=lambda x: (_ for _ in ()).throw(ValueError(f'Invalid number {x}')))
    if config.get('schema_version') != SCHEMA_VERSION:
        raise ValueError(f'Unsupported config schema_version {config.get("schema_version")!r} (this code reads '
                         f'{SCHEMA_VERSION}). A run archived by an older version replays with that version\'s code.')
    _keys(config, set(GROUPS) | SECTIONS | OPTIONAL_SECTIONS, 'config', set(GROUPS) | SECTIONS)
    # Backward-compatible default, materialized in every new run snapshot.
    config['simulation'].setdefault('nonlinear_max_evaluations', 20)
    def finite(value):
        if isinstance(value, dict):
            for v in value.values(): finite(v)
        elif isinstance(value, list):
            for v in value: finite(v)
        elif isinstance(value, float) and not math.isfinite(value):
            raise ValueError('All config numbers must be finite.')
    finite(config)
    values = {}
    for group, fields in GROUPS.items():
        extra = GROUP_EXTRAS.get(group, set())
        _keys(config[group], set(fields) | extra, group, set(fields) | extra)
        for key, (symbol, scale) in fields.items():
            values[symbol] = _scaled(config[group][key], scale)
    for key in GROUPS['files']:
        symbol = GROUPS['files'][key][0]
        values[symbol] = str((path.parent / config['files'][key]).resolve())
    for key in ('PROP_DT','PROP_SAVE_DT','PROP_RTOL','BANG_BANG_SAMPLE_PERIOD','MAX_THRUST_HOLD_TIME','COPV_HEADER_VOLUME','THROAT_AREA','EXIT_AREA','RAMP_BREAKPOINT_SPACING'):
        if values[key] <= 0:
            raise ValueError(f'{key} must be positive.')
    for key in ('MAIN_FUEL_OPEN_TIME','MAIN_OX_OPEN_TIME','OX_CLOSE_DELAY','FUEL_CLOSE_DELAY',
                'THROTTLE_STEP_TIME','TAIL_OFF_TIME'):
        if values[key] < 0:
            raise ValueError(f'{key} cannot be negative.')
    if values['MAIN_VALVE_CLOSE_TIME'] <= 0:
        raise ValueError('sequence.main_valve_close_s must be positive.')
    if max(values['MAIN_FUEL_OPEN_TIME'], values['MAIN_OX_OPEN_TIME']) >= values['PROFILE_START_TIME']:
        raise ValueError('sequence: OTV and FTV must start opening before main_ignition_s.')
    if values['THROTTLE_RATE_LIMIT'] is not None and values['THROTTLE_RATE_LIMIT'] <= 0:
        raise ValueError('throttle_rate_limit_lbf_s must be positive, or null for no limit.')
    for key in ('REGEN_PROPERTY_TEMPERATURE_POINTS','CHAMBER_MAP_PRESSURE_POINTS','CHAMBER_MAP_MIXTURE_RATIO_POINTS'):
        if not isinstance(values[key],int) or isinstance(values[key],bool) or values[key] < 2:
            raise ValueError(f'{key} must be an integer of at least two.')
    if not isinstance(values['PROP_MAX_NFEV'], int) or isinstance(values['PROP_MAX_NFEV'], bool) or values['PROP_MAX_NFEV'] < 1:
        raise ValueError('nonlinear_max_evaluations must be a positive integer.')
    for low, high in [('CHAMBER_MAP_PRESSURE_LOWER_LIMIT','CHAMBER_MAP_PRESSURE_UPPER_LIMIT'),
                      ('CHAMBER_MAP_MIXTURE_RATIO_LOWER_LIMIT','CHAMBER_MAP_MIXTURE_RATIO_UPPER_LIMIT')]:
        if not 0 < values[low] < values[high]:
            raise ValueError(f'{low} must be positive and less than {high}.')
    if values['COPV_MODEL'] not in ('lumped','header'):
        raise ValueError('copv_model must be lumped or header.')
    if values['PROP_PLOT_UNITS'] not in ('E','SI'):
        raise ValueError('plot_units must be E or SI.')
    if values['PROP_JACOBIAN_METHOD'] not in ('2-point','3-point'):
        raise ValueError('jacobian_method must be 2-point or 3-point.')
    for key,symbol in [('controller_steps_per_period','CONTROLLER_STEPS_PER_PERIOD'),('valve_steps_per_stroke','VALVE_STEPS_PER_STROKE')]:
        value = config['numerics'][key]
        if not isinstance(value,int) or isinstance(value,bool) or value < 1:
            raise ValueError(f'{key} must be a positive integer.')
        values[symbol] = value
    profile = config['thrust_profile']
    curves = _thrust_profiles(profile, path.parent)
    ramp = profile['ramp']
    values.update(PROP_THRUST_CURVE_NAME=profile['selected'], PROP_PROFILE_INTERPOLATION=profile['interpolation'],
                  PROP_THRUST_CURVE=curves[profile['selected']], THRUST_CURVES=curves,
                  THRUST_CURVE_MIN_SPACING=profile['csv_min_spacing_s'],
                  THRUST_CURVE_1DOF=curves.get('1dof'), RAMP_START_THRUST=ramp['start_lbf']*LBF,
                  RAMP_END_THRUST=ramp['end_lbf']*LBF, RAMP_DURATION=ramp['duration_s'], RAMP_POINTS=ramp['points'])
    for prefix,key in [('CSTAR','cstar_efficiency'),('CF','cf_efficiency')]:
        function, csv = efficiency(config['engine'][key], path.parent)
        values[prefix+'_EFFICIENCY_FUNCTION'], values[prefix+'_EFFICIENCY_CSV'] = function, csv
        values[prefix+'_EFFICIENCY_SCALE'] = 1.0   # the Monte Carlo moves the whole curve with this
    _keys(config['tanks'], ('oxidizer','fuel'), 'tanks', ('oxidizer','fuel'))
    for key,symbol in [('oxidizer','LOX_TANK'),('fuel','FUEL_TANK')]:
        tank = build_component(config['tanks'][key])
        if not isinstance(tank, elements.Tank) or tank.propellant != config['fluids'][key]:
            raise ValueError(f'{key} tank must be type Tank and match the configured fluid.')
        values[symbol] = tank
    values['COPVS'] = [build_component(b) for b in config['bottles']]
    if not values['COPVS'] or any(not isinstance(b,elements.COPVSource) for b in values['COPVS']):
        raise ValueError('At least one COPVSource is required.')
    paths = {'oxidizer_pressurant':'LOX_PRESSURANT_LEG','fuel_pressurant':'FUEL_PRESSURANT_LEG','oxidizer_feed':'LOX_FEED','fuel_feed':'FUEL_FEED'}
    _keys(config['flow_paths'], paths, 'flow_paths', paths)
    for key,symbol in paths.items():
        values[symbol] = [build_component(e) for e in config['flow_paths'][key]]
        if not values[symbol]:
            raise ValueError(f'{key} cannot be empty.')
        allowed = (elements.GasLine, elements.Orifice) if 'pressurant' in key else (elements.Line,elements.Orifice,elements.ThrottleValve,elements.CavitatingVenturi)
        if any(not isinstance(e,allowed) for e in values[symbol]):
            raise ValueError(f'{key}: unsupported component in this flow path.')

    # Side paths. schedule "torch" in a valve uses the firing window of the torch section.
    torch = config['torch']
    _keys(torch, ('enabled','open_time_s','duration_s','valve_stroke_s','branches','network','sequence','standalone_purge'),
          'torch', ('enabled','open_time_s','duration_s','valve_stroke_s'))
    if ('branches' in torch) == ('network' in torch):
        raise ValueError('torch: specify exactly one of network or legacy branches.')
    start, duration, stroke = torch['open_time_s'], torch['duration_s'], torch['valve_stroke_s']
    if (not isinstance(torch['enabled'], bool) or any(isinstance(x, bool) or not isinstance(x, (int,float))
            or not math.isfinite(x) for x in (start,duration,stroke))
            or start < 0 or duration <= 0 or stroke <= 0 or 2*stroke > duration):
        raise ValueError('torch: enabled must be boolean, start nonnegative, and duration >= 2 positive valve strokes.')
    torch_schedule = ([0.0, start, start + stroke, start + duration - stroke, start + duration], [0.0, 0.0, 1.0, 1.0, 0.0])
    schedules = {'torch': torch_schedule}
    values['REGULATORS'] = [build_component(r) for r in config['regulators']]
    if any(not isinstance(r, elements.Regulator) for r in values['REGULATORS']):
        raise ValueError('regulators: every entry must be type Regulator.')
    values['BRANCHES'] = [build_branch(b, schedules) for b in config['side_paths']]
    sources = {'ambient','Main Chamber','COPV','LOX Ullage','Fuel Ullage','LOX Tank Liquid','Fuel Tank Liquid'}
    sources.update(e.outlet_node for key in ('LOX_FEED','FUEL_FEED') for e in values[key][:-1])
    sources.update(r.outlet_node for r in values['REGULATORS'])
    values['TORCH_NETWORK'] = (validate_torch_network(torch['network'], sources) if 'network' in torch else None)
    values['TORCH_LEGS'], values['TORCH_SEQUENCE'] = None, None
    if values['TORCH_NETWORK'] is not None and torch['enabled']:
        values['TORCH_LEGS'] = _torch_legs(values['TORCH_NETWORK'])
        values['TORCH_SEQUENCE'] = _torch_sequence(torch.get('sequence'), values)
    values['TORCH_BRANCHES'] = [build_branch(b, schedules) for b in torch.get('branches', [])]
    values['NOT_MODELLED'] = dict(config['not_modelled'])
    values.update(TORCH_ENABLED=torch['enabled'], TORCH_OPEN_TIME=start, TORCH_DURATION=duration,
                  TORCH_VALVE_STROKE=stroke, TORCH_SCHEDULE=torch_schedule)
    if values['TORCH_LEGS'] is not None:
        values['NOT_MODELLED']['Torch spark plug and Pc check'] = ('Spark and the 100 ms torch Pc go/no-go check are not '
                                                                    'modelled; torch propellant lights on arrival')
        values['NOT_MODELLED']['OTHX heating'] = 'OTHX is a plain tube: no heating or vaporization of the torch LOX'
    elif torch['enabled']:
        values['BRANCHES'] = values['BRANCHES'] + values['TORCH_BRANCHES']
    elif values['TORCH_NETWORK'] is not None:
        values['NOT_MODELLED']['Torch fluids network'] = 'Configured fuel, oxidizer, OTHX, purge tee, chamber and throat; torch.enabled is false.'
    else:
        values['NOT_MODELLED']['OIMV / FIMV / HX-1'] = 'Torch igniter feeds, switched off (torch.enabled = false)'
    values.update(_purge(config['purge'], values['BRANCHES']))
    values.update(_monte_carlo(config['monte_carlo'], path.parent))

    objects = (values['COPVS'] + [values['LOX_TANK'],values['FUEL_TANK']] + sum([values[s] for s in paths.values()],[])
               + values['REGULATORS']
               + [e for b in values['BRANCHES'] + values['TORCH_BRANCHES'][:0 if torch['enabled'] and values['TORCH_LEGS'] is None else None] for e in b.elements])
    names = [e.name for e in objects]
    if values['TORCH_NETWORK'] is not None:
        network = values['TORCH_NETWORK']
        names += [network['chamber']['name'], network['throat']['name']]
        names += [e['name'] for key in ('fuel_leg','purge_leg','oxidizer_leg') for e in network[key]['elements']]
    if len(set(names)) != len(names):
        raise ValueError(f'All source and component names must be unique: {sorted({n for n in names if names.count(n) > 1})}.')
    values.update(EXPANSION_RATIO=values['EXIT_AREA']/values['THROAT_AREA'],
                  COPV_COUNT=len(values['COPVS']), COPV_VOLUME=sum(b.volume for b in values['COPVS']),
                  COPV_UNIT_VOLUME=values['COPVS'][0].volume,
                  PRESSURANT_INITIAL_PRESSURE=values['COPVS'][0].pressure,
                  PRESSURANT_INITIAL_TEMPERATURE=values['COPVS'][0].temperature,
                  COPV_MIN_TEMPERATURES=[b.min_temperature for b in values['COPVS']])
    return config, values, path


def integration_step(values, requested=None):
    requested = values['PROP_DT'] if requested is None else requested
    period = values['BANG_BANG_SAMPLE_PERIOD']
    travel = [t for key in ('LOX_FEED','FUEL_FEED','LOX_PRESSURANT_LEG','FUEL_PRESSURANT_LEG')
              for e in values[key] for t in (getattr(e,'open_time',0),getattr(e,'close_time',0)) if t > 0]
    step = min(requested, period/values['CONTROLLER_STEPS_PER_PERIOD'],
               min(travel)/values['VALVE_STEPS_PER_STROKE'] if travel else requested)
    return period / math.ceil(period/step)
