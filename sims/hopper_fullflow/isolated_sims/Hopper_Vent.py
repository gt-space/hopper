"""Tank vent time: how long each vent Branch in HOPPER_DATA.py takes to blow its
tank down, as a small FullFlow transient per vent.

The ullage is a FullFlow Volume of pressurant gas, vented through the branch's
restrictions (one SeriesCdA -> CompressibleOrifice to ambient). Two cases per
tank, both starting at the tank set pressure (+ upper limit):
    full    ullage = tank volume - loaded liquid (pad pressure before the burn)
    empty   ullage = whole tank volume (after the burn)
The gas expands adiabatically (no heat from the walls or liquid) and LOX
boil-off during the vent is not included; VENT_* inputs are in HOPPER_DATA.py.

    python isolated_sims/Hopper_Vent.py      prints the table and writes vent_times.csv
"""

import csv
import os

import h5py
import numpy as np
from thermoprop import Fluid
from fullflow import CompressibleOrifice, Lookup, Network, SeriesCdA, State, Transient, Volume
from scipy.constants import liter, psi

import sys
from pathlib import Path

if not __package__:   # run as a script (python isolated_sims/Hopper_Vent.py): make HOPPER_DATA and the Backend
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))   # package importable

import HOPPER_DATA as D
from Backend.Hopper_Prop_System import liquid_properties


def vent_time(branch, tank, ullage_volume, filename):
    """Time [s] for the ullage to fall from tank set + upper limit to
    AMBIENT + VENT_END_GAUGE_PRESSURE, and the pressure history."""
    start = tank.set_pressure + tank.upper_limit
    end = D.AMBIENT_PRESSURE + D.VENT_END_GAUGE_PRESSURE
    network = Network("Vent")
    gas = Lookup("Ullage Gas", network, Fluid, D.PRESSURANT, pressure=start,
                 temperature=tank.ullage_temperature, priority=("enthalpy", "temperature"))
    path = SeriesCdA("Vent Series CdA", network, effective_areas=[e.cda_value for e in branch.elements])
    orifice = CompressibleOrifice("Vent", network, upstream_total_pressure=gas.pressure,
                                  upstream_total_temperature=gas.temperature, downstream_pressure=D.AMBIENT_PRESSURE,
                                  discharge_coefficient=1.0, cross_sectional_area=path.effective_area,
                                  gas_constant=gas.gas_constant, specific_heat_ratio=gas.gamma)
    Volume("Ullage", network, pressure=gas.pressure, temperature=gas.temperature, volume=ullage_volume,
           density=gas.density, internal_energy=gas.internal_energy, enthalpy=gas.enthalpy,
           energy_variable="enthalpy", mass_flow_out=State.maximum(orifice.mass_flow, 0.0))
    network.track("Pressure [psia]", gas.pressure / psi)
    network.track("Temperature [K]", gas.temperature)
    Transient(network).solve(dt=D.VENT_DT, t_final=D.VENT_MAX_TIME, save_dt=D.VENT_DT, filename=filename,
                             rtol=D.PROP_RTOL, verbose=False, statistics=False)
    with h5py.File(filename, "r") as h5:
        run = h5["Vent/transient/runs/base"]
        time = np.asarray(run["time"][()], float)
        pressure = np.asarray(run["tracks/Pressure_[psia]"][()], float)
        temperature = np.asarray(run["tracks/Temperature_[K]"][()], float)
    below = np.nonzero(pressure * psi <= end)[0]
    if below.size == 0:
        return float("nan"), time, pressure, temperature
    k = max(below[0], 1)   # linear crossing between the saved points either side
    t_end = float(time[k - 1] + (end / psi - pressure[k - 1]) * (time[k] - time[k - 1]) / (pressure[k] - pressure[k - 1]))
    return t_end, time, pressure, temperature


def vent_times(directory=D.PROP_RESULTS_DIRECTORY):
    """One row per vent branch and case; writes vent_times.csv to directory."""
    tanks = {"LOX Ullage": D.LOX_TANK, "Fuel Ullage": D.FUEL_TANK}
    os.makedirs(directory, exist_ok=True)
    scratch = os.path.join(directory, "_vent.h5")
    rows = []
    for branch in D.BRANCHES:
        tank = tanks.get(branch.source)
        if tank is None or not any(getattr(e, "schedule", None) is not None for e in branch.elements):
            continue   # only commanded vents off a tank ullage (reliefs are not vents)
        density = liquid_properties(tank.propellant, tank.set_pressure, tank.liquid_temperature)["density"]
        cases = {"full": tank.volume - tank.liquid_mass / density, "empty": tank.volume}
        for case, volume in cases.items():
            t_end, *_ = vent_time(branch, tank, volume, scratch)
            requirement = D.VENT_TIME_REQUIREMENTS.get(branch.name)
            rows.append({"Vent": branch.name, "Tank": tank.name, "Case": case,
                         "Gas volume [L]": round(volume / liter, 2),
                         "CdA [mm2]": round(1e6 / sum(1 / e.cda_value**2 for e in branch.elements) ** 0.5, 2),
                         "From [psia]": round((tank.set_pressure + tank.upper_limit) / psi, 1),
                         "To [psia]": round((D.AMBIENT_PRESSURE + D.VENT_END_GAUGE_PRESSURE) / psi, 1),
                         "Vent time [s]": round(t_end, 1),
                         "Requirement [s]": requirement if requirement is not None else "",
                         "Meets": "" if requirement is None else ("yes" if t_end <= requirement else "NO")})
    if os.path.exists(scratch):
        os.remove(scratch)
    path = os.path.join(directory, "vent_times.csv")
    with open(path, "w", newline="", encoding="utf-8") as file:
        writer = csv.DictWriter(file, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return rows, path


if __name__ == "__main__":
    rows, path = vent_times()
    for row in rows:
        print(row)
    print("wrote", path)
