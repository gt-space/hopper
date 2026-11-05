"""JSON configuration boundary. Physics continues to consume SI-valued inputs."""

import inspect
import json
import math
import os
from pathlib import Path

from scipy.constants import psi, inch
import Hopper_Elements as elements
from Hopper_Inputs import ramp_profile, validate_profile, table_value

# Unit conversions and solver bindings are code, not vehicle settings.
from Hopper_Config_Fields import GROUPS, UNITS

ROOT = Path(__file__).resolve().parent
TYPES = {name: getattr(elements, name) for name in (
    'Tank', 'COPVSource', 'Line', 'GasLine', 'Orifice', 'Injector', 'Regen',
    'Valve', 'BangBangValve', 'CheckValve', 'ThrottleValve', 'CavitatingVenturi')}


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


def build_component(spec):
    kind = spec.get('type')
    if kind not in TYPES:
        raise ValueError(f'Unknown component type {kind!r}; choose from {list(TYPES)}.')
    cls = TYPES[kind]
    parameters = inspect.signature(cls.__init__).parameters
    allowed = {'type'}
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
        elif value is not None and scale != 1:
            value *= scale
        kwargs[name] = value
    _keys(spec, allowed, spec.get('name', kind))
    for field in ('open_time', 'close_time', 'roughness', 'ullage_htc'):
        if field in kwargs and kwargs[field] < 0:
            raise ValueError(f'{kwargs["name"]}: {field} cannot be negative.')
    for field in ('volume', 'pressure', 'liquid_mass', 'liquid_temperature', 'ullage_temperature',
                  'radius', 'set_pressure', 'cda', 'max_cda', 'throat_area'):
        if field in kwargs and kwargs[field] <= 0:
            raise ValueError(f'{kwargs["name"]}: {field} must be positive.')
    for field in ('upper_limit', 'lower_limit'):
        if field in kwargs and kwargs[field] < 0:
            raise ValueError(f'{kwargs["name"]}: deadband cannot be negative.')
    schedule = kwargs.get('schedule')
    if schedule is not None:
        times, openings = schedule
        if not times or len(times) != len(openings) or any(b <= a for a,b in zip(times,times[1:])) or any(not 0 <= x <= 1 for x in openings):
            raise ValueError(f'{kwargs["name"]}: invalid valve schedule [times_s, openings].')
    return cls(**kwargs)


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


def load_config(filename=None):
    path = Path(filename or os.environ.get('HOPPER_CONFIG', ROOT / 'hopper_config.json')).resolve()
    config = json.loads(path.read_text(encoding='utf-8-sig'), object_pairs_hook=_unique_object,
                        parse_constant=lambda x: (_ for _ in ()).throw(ValueError(f'Invalid number {x}')))
    _keys(config, set(GROUPS) | {'schema_version','description','thrust_profile','bottles','tanks','flow_paths'}, 'config',
          set(GROUPS) | {'schema_version','thrust_profile','bottles','tanks','flow_paths'})
    if config['schema_version'] != 1:
        raise ValueError('Unsupported config schema_version.')
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
        extra = {'engine':{'cstar_efficiency','cf_efficiency'}, 'numerics':{'controller_steps_per_period','valve_steps_per_stroke'}}.get(group,set())
        _keys(config[group], set(fields) | extra, group, set(fields) | extra)
        for key, (symbol, scale) in fields.items():
            raw = config[group][key]
            values[symbol] = raw * scale if scale != 1 else raw
    for key in ('chamber_map','history','results_directory'):
        symbol = GROUPS['files'][key][0]
        values[symbol] = str((path.parent / config['files'][key]).resolve())
    for key in ('PROP_DT','PROP_SAVE_DT','PROP_RTOL','BANG_BANG_SAMPLE_PERIOD','MAX_THRUST_HOLD_TIME','COPV_HEADER_VOLUME','THROAT_AREA','EXIT_AREA'):
        if values[key] <= 0:
            raise ValueError(f'{key} must be positive.')
    for key in ('MAIN_VALVE_OPEN_TIME','MAIN_VALVE_RAMP_TIME','THROTTLE_STEP_TIME','TAIL_OFF_TIME'):
        if values[key] < 0:
            raise ValueError(f'{key} cannot be negative.')
    for key in ('REGEN_PROPERTY_TEMPERATURE_POINTS','CHAMBER_MAP_PRESSURE_POINTS','CHAMBER_MAP_MIXTURE_RATIO_POINTS'):
        if not isinstance(values[key],int) or isinstance(values[key],bool) or values[key] < 2:
            raise ValueError(f'{key} must be an integer of at least two.')
    if not isinstance(values['PROP_MAX_NFEV'], int) or isinstance(values['PROP_MAX_NFEV'], bool) or values['PROP_MAX_NFEV'] < 1:
        raise ValueError('nonlinear_max_evaluations must be a positive integer.')
    for low, high in [('CHAMBER_MAP_PRESSURE_LOWER_LIMIT','CHAMBER_MAP_PRESSURE_UPPER_LIMIT'),
                      ('CHAMBER_MAP_MIXTURE_RATIO_LOWER_LIMIT','CHAMBER_MAP_MIXTURE_RATIO_UPPER_LIMIT')]:
        if not 0 < values[low] < values[high]:
            raise ValueError(f'{low} must be positive and less than {high}.')
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
    _keys(profile, ('selected','interpolation','ramp','profiles'), 'thrust_profile', ('selected','interpolation','ramp','profiles'))
    if profile['interpolation'] not in ('linear','step'):
        raise ValueError('Thrust interpolation must be linear or step.')
    ramp = profile['ramp']
    _keys(ramp, ('start_lbf','end_lbf','duration_s','points'), 'ramp', ('start_lbf','end_lbf','duration_s','points'))
    curves = {}
    for name, entry in profile['profiles'].items():
        if name == 'ramp':
            raise ValueError('The profile name ramp is reserved for the configured ramp.')
        _keys(entry, ('time_s_thrust_n',), name, ('time_s_thrust_n',))
        pairs = entry['time_s_thrust_n']
        if len(dict(pairs)) != len(pairs):
            raise ValueError(f'{name}: duplicate thrust times.')
        curves[name] = validate_profile(dict(pairs))
    from Hopper_Inputs import LBF
    curves['ramp'] = ramp_profile(ramp['start_lbf']*LBF, ramp['end_lbf']*LBF, ramp['duration_s'], ramp['points'])
    if profile['selected'] not in curves:
        raise ValueError('Selected thrust profile is not defined.')
    values.update(PROP_THRUST_CURVE_NAME=profile['selected'], PROP_PROFILE_INTERPOLATION=profile['interpolation'],
                  PROP_THRUST_CURVE=curves[profile['selected']], THRUST_CURVES=curves,
                  THRUST_CURVE_1DOF=curves.get('1dof'), RAMP_START_THRUST=ramp['start_lbf']*LBF,
                  RAMP_END_THRUST=ramp['end_lbf']*LBF, RAMP_DURATION=ramp['duration_s'], RAMP_POINTS=ramp['points'])
    for prefix,key in [('CSTAR','cstar_efficiency'),('CF','cf_efficiency')]:
        function, csv = efficiency(config['engine'][key], path.parent)
        values[prefix+'_EFFICIENCY_FUNCTION'], values[prefix+'_EFFICIENCY_CSV'] = function, csv
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
    objects = values['COPVS'] + [values['LOX_TANK'],values['FUEL_TANK']] + sum([values[s] for s in paths.values()],[])
    names = [e.name for e in objects]
    if len(set(names)) != len(names):
        raise ValueError('All source and component names must be unique.')
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
