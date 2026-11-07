"""Archive each run's input JSON, referenced assets and resolved solver settings."""

import copy
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import shutil


# Inputs that stay in the archive when the results are copied to the viewing folder.
_ARCHIVE_ONLY = {'input_config.json', 'run_config.json', 'resolved_inputs.json'}


def prepare_run(data):
    """Call before importing physics modules so their bindings use archived inputs."""
    from Backend.Hopper_Config import integration_step
    root = Path(data.PROP_ARCHIVE_DIRECTORY)
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
    for name, entry in config['thrust_profile']['profiles'].items():
        if 'csv' in entry:
            _, entry['csv'] = archive(Path(data.CONFIG_PATH).parent / entry['csv'], f'thrust_{name}.csv')
    map_path = assets / 'chamber_map.h5'
    if not data.PROP_GENERATE_CHAMBER_MAP and Path(data.CHAMBER_MAP_FILENAME).exists():
        archive(data.CHAMBER_MAP_FILENAME, map_path.name)
    data.CHAMBER_MAP_FILENAME = str(map_path.resolve())
    # finish_run copies the results from the archive to where they are viewed.
    data.PROP_LATEST_DIRECTORY = data.PROP_RESULTS_DIRECTORY
    data.PROP_LATEST_FILENAME = data.PROP_FILENAME
    data.PROP_FILENAME = str((output / Path(data.PROP_FILENAME).name).resolve())
    data.PROP_RESULTS_DIRECTORY = str(output.resolve())
    config['files'] = {'chamber_map':'inputs/chamber_map.h5', 'history':Path(data.PROP_FILENAME).name,
                       'results_directory':'.', 'run_archive_directory':'.'}
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


def finish_run(data):
    """Copy a finished run's reports and plots to the viewing folder (Outputs), and move
    its raw history there, so one large history is kept rather than one per run.
    outputs.keep_run_history = true also leaves a copy in the archive."""
    archive, latest = Path(data.PROP_RESULTS_DIRECTORY).resolve(), Path(data.PROP_LATEST_DIRECTORY).resolve()
    if archive == latest:
        return
    latest.mkdir(parents=True, exist_ok=True)
    for file in archive.iterdir():
        if file.is_file() and file.suffix != '.h5' and file.name not in _ARCHIVE_ONLY:
            shutil.copy2(file, latest / file.name)
    history, target = Path(data.PROP_FILENAME), Path(data.PROP_LATEST_FILENAME)
    if not history.exists():
        return
    try:
        target.parent.mkdir(parents=True, exist_ok=True)
        if data.PROP_KEEP_RUN_HISTORY:
            shutil.copy2(history, target)
        else:
            os.replace(history, target)
    except PermissionError:   # e.g. the old history is open in another program
        print(f'Could not replace {target}; the raw history stays at {history}.')
    print(f'Latest results copied to {latest}')
