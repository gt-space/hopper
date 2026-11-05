"""Archive each run's input JSON, referenced assets and resolved solver settings."""

import copy
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import shutil


def prepare_run(data):
    """Call before importing physics modules so their bindings use archived inputs."""
    from Hopper_Config import integration_step
    root = Path(data.PROP_RESULTS_DIRECTORY)
    output = root / datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S_%fZ')
    output.mkdir(parents=True, exist_ok=False)
    config = copy.deepcopy(data.CONFIG)
    (output / 'input_config.json').write_text(data.CONFIG_SOURCE_TEXT, encoding='utf-8')
    assets = output / 'inputs'
    assets.mkdir()
    hashes = {}

    def archive(filename, name):
        source = Path(filename)
        target = assets / name
        shutil.copy2(source, target)
        hashes[str(target.relative_to(output))] = hashlib.sha256(target.read_bytes()).hexdigest()
        return str(target.resolve()), target.relative_to(output).as_posix()

    for prefix, key in [('CSTAR','cstar_efficiency'),('CF','cf_efficiency')]:
        filename = getattr(data, prefix + '_EFFICIENCY_CSV')
        if filename:
            absolute, relative = archive(filename, key + '.csv')
            setattr(data, prefix + '_EFFICIENCY_CSV', absolute)
            config['engine'][key]['csv'] = relative
    map_path = assets / 'chamber_map.h5'
    if not data.PROP_GENERATE_CHAMBER_MAP and Path(data.CHAMBER_MAP_FILENAME).exists():
        archive(data.CHAMBER_MAP_FILENAME, map_path.name)
    data.CHAMBER_MAP_FILENAME = str(map_path.resolve())
    data.PROP_FILENAME = str((output / Path(data.PROP_FILENAME).name).resolve())
    data.PROP_RESULTS_DIRECTORY = str(output.resolve())
    config['files'] = {'chamber_map':'inputs/chamber_map.h5', 'history':Path(data.PROP_FILENAME).name, 'results_directory':'.'}
    data.RUN_CONFIG = config
    (output / 'run_config.json').write_text(json.dumps(config,indent=2)+'\n', encoding='utf-8')
    resolved = {'source_config':str(data.CONFIG_PATH), 'config':config,
        'derived':{'throat_diameter_in':(4*data.THROAT_AREA/3.141592653589793)**.5/.0254,
                   'expansion_ratio':data.EXPANSION_RATIO,
                   'integration_step_cap_s':integration_step(vars(data)),
                   'selected_profile_time_s_thrust_n':list(data.PROP_THRUST_CURVE.items()),
                   'header_initial_pressure_pa':data.PRESSURANT_INITIAL_PRESSURE,
                   'header_initial_temperature_k':data.PRESSURANT_INITIAL_TEMPERATURE},
        'assets_sha256':hashes}
    data.RESOLVED_INPUTS = resolved
    write_resolved(data)
    print(f'Run inputs and results: {output}')
    return output


def write_resolved(data):
    """Refresh map provenance after optional generation; called before the solve."""
    resolved = data.RESOLVED_INPUTS
    path = Path(data.CHAMBER_MAP_FILENAME)
    if path.exists():
        resolved['assets_sha256']['inputs/chamber_map.h5'] = hashlib.sha256(path.read_bytes()).hexdigest()
    Path(data.PROP_RESULTS_DIRECTORY, 'resolved_inputs.json').write_text(
        json.dumps(resolved,indent=2)+'\n', encoding='utf-8')
