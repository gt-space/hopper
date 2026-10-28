"""LOX drain calculator: how long the loaded LOX takes to drain through OFD.

One small FullFlow transient per ullage pressure in LOX_DRAIN_ULLAGE_PRESSURES:
a DrainingTank (Hopper_Components.py: liquid mass state, outlet pressure =
ullage + liquid head) feeding the drain path, i.e. the LOX_DRAIN_BRANCH
restrictions plus LOX_DRAIN_EXTRA_ELEMENTS as one FullFlow SeriesCdA into a
DischargeCoefficient to ambient. The ullage is held at its pressure for the
whole drain (vented, or a regulated GN2 push). OFD is the same valve as OVNT.

Assumptions: incompressible LOX at the loaded temperature, no boil-off, no
flashing in the valve. A vented (0 psig) drain of near-saturated LOX will
partly flash, so its time is optimistic; the subcooling column is how far the
outlet pressure sits above the LOX vapour pressure at the start.

    python isolated_sims/Hopper_Drain.py     prints the table and writes lox_drain_times.csv
"""

import csv
import math
import os

import h5py
import numpy as np
from fullflow import DischargeCoefficient, Network, SeriesCdA, State, Transient
from scipy.constants import g, psi

import sys
from pathlib import Path

if not __package__:   # run as a script (python isolated_sims/Hopper_Drain.py): make HOPPER_DATA and the Backend
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))   # package importable

import HOPPER_DATA as D
from Backend.Hopper_Components import DrainingTank
from Backend.Hopper_Prop_System import liquid_properties


def drain_path():
    """Restrictions in series from the tank bottom to the drain exit."""
    branch = next((b for b in D.BRANCHES if b.name == D.LOX_DRAIN_BRANCH), None)
    if branch is None:
        raise ValueError(f"LOX_DRAIN_BRANCH {D.LOX_DRAIN_BRANCH!r} is not in BRANCHES.")
    return list(branch.elements) + list(D.LOX_DRAIN_EXTRA_ELEMENTS)


def drain_cda():
    """Series CdA of the drain path [m^2]."""
    return sum(1.0 / element.cda_value**2 for element in drain_path()) ** -0.5


def drain_time(ullage_gauge_pressure, filename, liquid_mass=None):
    """(drain time [s], time, liquid mass, mass flow) for one ullage gauge pressure."""
    tank = D.LOX_TANK
    liquid_mass = tank.liquid_mass if liquid_mass is None else liquid_mass
    ullage = D.AMBIENT_PRESSURE + ullage_gauge_pressure
    density = liquid_properties(tank.propellant, ullage, tank.liquid_temperature)["density"]
    area = tank.volume / D.LOX_TANK_HEIGHT
    empty = D.LOX_DRAIN_END_FRACTION * liquid_mass

    # Run length: twice the time at the starting flow covers a pure gravity drain
    # (whose flow falls to zero as the head goes) with margin.
    head = liquid_mass / (density * area) + D.LOX_DRAIN_OUTLET_DROP
    start_flow = drain_cda() * math.sqrt(2.0 * density * (ullage_gauge_pressure + density * g * head))
    t_final = min(D.DRAIN_MAX_TIME, 2.5 * liquid_mass / start_flow)

    network = Network("LOX Drain")
    mass = State(liquid_mass)
    tank_node = DrainingTank("LOX Tank", network, liquid_mass=mass, ullage_pressure=ullage, density=density,
                             cross_sectional_area=area, outlet_drop=D.LOX_DRAIN_OUTLET_DROP)
    path = SeriesCdA("Drain Series CdA", network, effective_areas=[e.cda_value for e in drain_path()])
    # The opening fades out over the last LOX_DRAIN_END_FRACTION of the load so the run settles at empty.
    opening = (mass / empty).clip(0.0, 1.0)
    valve = DischargeCoefficient("Drain", network, upstream_pressure=tank_node.outlet_pressure,
                                 downstream_pressure=D.AMBIENT_PRESSURE, density=density,
                                 discharge_coefficient=opening, cross_sectional_area=path.effective_area)
    tank_node.mass_flow_out = State.maximum(valve.mass_flow, 0.0)
    network.track("Liquid Mass [kg]", mass)
    network.track("Drain Mass Flow [kg/s]", tank_node.mass_flow_out)
    dt = t_final / D.DRAIN_STEPS
    Transient(network).solve(dt=dt, t_final=t_final, save_dt=dt, filename=filename, rtol=D.PROP_RTOL,
                             verbose=False, statistics=False)

    with h5py.File(filename, "r") as h5:
        run = h5["LOX_Drain/transient/runs/base"]
        time = np.asarray(run["time"][()], float)
        masses = np.asarray(run["tracks/Liquid_Mass_[kg]"][()], float)
        flows = np.asarray(run["tracks/Drain_Mass_Flow_[kg_s]"][()], float)
    target = 1.01 * empty   # the fade-out approaches "empty" asymptotically
    below = np.nonzero(masses <= target)[0]
    if below.size == 0:
        return float("nan"), time, masses, flows
    k = max(below[0], 1)
    t_end = float(time[k - 1] + (target - masses[k - 1]) * (time[k] - time[k - 1]) / (masses[k] - masses[k - 1]))
    return t_end, time, masses, flows


def drain_times(directory=D.PROP_RESULTS_DIRECTORY):
    """One row per ullage pressure; writes lox_drain_times.csv to directory."""
    os.makedirs(directory, exist_ok=True)
    scratch = os.path.join(directory, "_drain.h5")
    tank = D.LOX_TANK
    requirement = D.LOX_DRAIN_TIME_REQUIREMENT
    rows = []
    for gauge in D.LOX_DRAIN_ULLAGE_PRESSURES:
        t_end, time, masses, flows = drain_time(gauge, scratch)
        properties = liquid_properties(tank.propellant, D.AMBIENT_PRESSURE + gauge, tank.liquid_temperature)
        density = properties["density"]
        head = tank.liquid_mass / (density * tank.volume / D.LOX_TANK_HEIGHT) + D.LOX_DRAIN_OUTLET_DROP
        outlet = D.AMBIENT_PRESSURE + gauge + density * g * head
        rows.append({
            "Ullage [psig]": round(gauge / psi, 1),
            "LOX [kg]": tank.liquid_mass,
            "Drain CdA [mm2]": round(drain_cda() * 1e6, 2),
            "Liquid height [m]": round(head - D.LOX_DRAIN_OUTLET_DROP, 3),
            "Start flow [kg/s]": round(float(flows[0]), 4),
            "Drain time [s]": round(t_end, 1),
            "Drain time [min]": round(t_end / 60.0, 2),
            "Start subcooling [psi]": round((outlet - properties["vapor_pressure"]) / psi, 2),
            "Requirement [s]": "" if requirement is None else requirement,
            "Meets": "" if requirement is None else ("yes" if t_end <= requirement else "NO"),
        })
    if os.path.exists(scratch):
        os.remove(scratch)
    path = os.path.join(directory, "lox_drain_times.csv")
    with open(path, "w", newline="", encoding="utf-8") as file:
        writer = csv.DictWriter(file, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return rows, path


if __name__ == "__main__":
    rows, path = drain_times()
    for row in rows:
        print(row)
    print("wrote", path)
