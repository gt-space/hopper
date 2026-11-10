"""Purge flow check: can the GN2 purges match the propellant volumetric flow?

From the last full run (PROP_FILENAME) it takes
  * each propellant's volumetric flow along the burn, Q = mdot / rho, and its
    peak;
  * the COPV state when the purges open, at shutdown (end of life, the weakest supply);
  * the pressure in each feed line where its purge enters (OTV Outlet for the
    LOX purge, FTV Outlet for the fuel purge, from the P&ID).
Then, for each purge line in side_paths ("LOX Purge", "Fuel Purge": tube ->
valve -> check valve), it computes the GN2 flow the line would deliver into
the feed line at every instant of the burn, and that flow's volume at the
feed-line pressure, so the two volumetric flows can be compared.

Regulator (regulators, PR-REG): ideal pressure-reducing regulator. Its outlet
is held at min(COPV, ambient + set) and its own flow capacity is its CdA,
choked from the COPV state; the purge gets the smaller of what its line passes
and what the regulator can supply. The gas is throttled at constant enthalpy
(real-gas N2, ThermoProp) through the regulator and the valves.
Gas restrictions are one SeriesCdA through FullFlow's CompressibleOrifice
equations (ideal-gas flow with the real-gas R and gamma), as in the transient.

Outputs in PROP_RESULTS_DIRECTORY:
    purge_flow.csv   time series: propellant and purge volumetric flow, feed pressure
    purge_flow.png   propellant vs purge volumetric flow along the burn, per side

    python isolated_sims/Hopper_Purge.py
"""

import csv
import math
import os
import sys
from pathlib import Path

import numpy as np

if not __package__:   # run as a script (python isolated_sims/Hopper_Purge.py): make HOPPER_DATA and the Backend
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))   # package importable

import HOPPER_DATA as D
from scipy.constants import liter, milli, psi
from thermoprop import Fluid

from Backend.Hopper_Export import LoadTracks, _track
from Backend.Hopper_Report import read_metadata

# (side, purge branch, feed track for the propellant flow, its density track, feed node the purge enters)
SIDES = (("LOX", "LOX Purge", "Ox Mass Flow [kg/s]", "LOX Liquid Density [kg/m3]", "OTV Outlet"),
         ("Fuel", "Fuel Purge", "Fuel Mass Flow [kg/s]", "Fuel Liquid Density [kg/m3]", "FTV Outlet"))


def orifice_mass_flow(cda, p_up, t_up, p_down, gamma, gas_constant):
    """Ideal-gas orifice flow [kg/s], choked or subcritical (FullFlow CompressibleOrifice's equations)."""
    if p_up <= p_down:
        return 0.0
    ratio = p_down / p_up
    critical = (2.0 / (gamma + 1.0)) ** (gamma / (gamma - 1.0))
    if ratio <= critical:
        flow = math.sqrt(gamma / (gas_constant * t_up) * (2.0 / (gamma + 1.0)) ** ((gamma + 1.0) / (gamma - 1.0)))
    else:
        flow = math.sqrt(2.0 * gamma / (gas_constant * t_up * (gamma - 1.0))
                         * max(ratio ** (2.0 / gamma) - ratio ** ((gamma + 1.0) / gamma), 0.0))
    return cda * p_up * flow


def purge_start_time(filename, time):
    """Simulation time [s] of the shutdown request, when the purge timing starts (end of the
    run for an old file without metadata)."""
    metadata = read_metadata(filename)
    if metadata is None:
        return float(time[-1])
    return metadata["profile_end"]


def series_cda(elements):
    return sum(1.0 / element.cda_value**2 for element in elements) ** -0.5


def purge_check(filename=D.PROP_FILENAME, directory=D.PROP_RESULTS_DIRECTORY):
    """Write purge_flow.csv / .png from a results file; returns {metric: value}."""
    import matplotlib.pyplot as plt

    time, tracks = LoadTracks(filename)
    regulator = D.REGULATORS[0]
    branches = {b.name: b for b in D.BRANCHES}

    # End-of-life COPV and the regulator outlet (isenthalpic throttle).
    start = purge_start_time(filename, time)
    p_copv = float(np.interp(start, time, _track(tracks, "COPV Pressure [psia]"))) * psi
    t_copv = float(np.interp(start, time, _track(tracks, "COPV Temperature [K]")))
    copv = Fluid(D.PRESSURANT, pressure=p_copv, temperature=t_copv)
    p_reg = min(p_copv, D.AMBIENT_PRESSURE + regulator.set_pressure)
    reg_gas = Fluid(D.PRESSURANT, pressure=p_reg, enthalpy=copv.enthalpy)
    reg_capacity = (orifice_mass_flow(regulator.cda_value, p_copv, t_copv, p_reg, copv.gamma, copv.gas_constant)
                    if regulator.cda_value else math.inf)

    burning = _track(tracks, "Total Mass Flow [kg/s]") > 0.05 * _track(tracks, "Total Mass Flow [kg/s]").max()
    summary = {"COPV end-of-life pressure [psia]": round(p_copv / psi, 1),
               "COPV end-of-life temperature [K]": round(t_copv, 1),
               "Regulator outlet [psia]": round(p_reg / psi, 1),
               "Regulator outlet temperature [K]": round(reg_gas.temperature, 1),
               "Regulator capacity at EOL [g/s]": round(1000 * reg_capacity, 2)}
    columns = {"Time [s]": time}
    results = {}
    for side, branch_name, flow_track, density_track, node in SIDES:
        branch = branches[branch_name]
        cda = series_cda(branch.elements)
        q_prop = _track(tracks, flow_track) / _track(tracks, density_track)                 # [m^3/s]
        p_feed = _track(tracks, f"{node} Pressure [Pa]")
        mdot = np.array([min(orifice_mass_flow(cda, p_reg, reg_gas.temperature, p, reg_gas.gamma, reg_gas.gas_constant),
                             reg_capacity) for p in p_feed])
        rho_gas = np.array([Fluid(D.PRESSURANT, pressure=max(p, 1000.0), enthalpy=reg_gas.enthalpy).density for p in p_feed])
        q_purge = mdot / rho_gas
        k_peak = int(np.argmax(np.where(burning, q_prop, -1.0)))
        ratio = np.where(burning & (q_prop > 0), q_purge / np.maximum(q_prop, 1e-12), np.nan)
        results[side] = (q_prop, q_purge, p_feed, k_peak)
        columns.update({f"{side} propellant flow [L/s]": q_prop / liter, f"{side} purge GN2 flow [L/s]": q_purge / liter,
                        f"{side} purge GN2 mass flow [g/s]": 1000 * mdot, f"{node} pressure [psia]": p_feed / psi})
        summary.update({
            f"{side} purge path CdA [mm2]": round(cda / milli**2, 2),
            f"{side} max propellant volumetric flow [L/s]": round(q_prop[k_peak] / liter, 3),
            f"{side} purge flow at that moment [L/s]": round(q_purge[k_peak] / liter, 3),
            f"{side} purge / propellant at peak [%]": round(100 * q_purge[k_peak] / q_prop[k_peak], 1),
            f"{side} feed pressure at peak [psia]": round(p_feed[k_peak] / psi, 1),
            f"{side} purge / propellant, worst point in burn [%]": round(100 * np.nanmin(ratio), 1),
            f"{side} purge meets max propellant flow": "yes" if q_purge[k_peak] >= q_prop[k_peak] else "NO",
        })
    both = sum(columns[f"{s} purge GN2 mass flow [g/s]"][results[s][3]] for s, *_ in SIDES) / 1000
    summary["Both purges at once vs regulator capacity [%]"] = (round(100 * both / reg_capacity, 1)
                                                                 if reg_capacity < math.inf else "no limit")

    os.makedirs(directory, exist_ok=True)
    names = list(columns)
    with open(os.path.join(directory, "purge_flow.csv"), "w", newline="", encoding="utf-8") as file:
        writer = csv.writer(file)
        writer.writerow(names)
        for i in range(len(time)):
            writer.writerow([round(float(columns[n][i]), 5) for n in names])

    fig, axes = plt.subplots(2, 1, figsize=(11, 8.5), sharex=True)
    for ax, (side, branch_name, *_ , node) in zip(axes, SIDES):
        q_prop, q_purge, p_feed, k = results[side]
        colour = "tab:blue" if side == "LOX" else "tab:red"
        ax.plot(time, q_prop / liter, color=colour, lw=1.5, label=f"{side} propellant flow")
        ax.plot(time, np.where(burning, q_purge, np.nan) / liter, color="k", lw=1.5,
                label=f"GN2 purge flow at {node} pressure")
        ax.axhline(q_prop[k] / liter, color=colour, ls="--", lw=1, label=f"peak propellant flow {q_prop[k] / liter:.3f} L/s")
        ax.plot(time[k], q_purge[k] / liter, "o", color="k")
        twin = ax.twinx(); twin.plot(time, p_feed / psi, color="0.75", lw=1); twin.set_ylabel(f"{node} pressure [psia]", color="0.5")
        ax.set_ylim(0, 2.0 * q_prop[k] / liter)   # start-up / shutdown (feed near ambient) would dwarf the burn
        ax.set_ylabel("Volumetric flow [L/s]"); ax.grid(alpha=0.3); ax.legend(fontsize=9, loc="lower center")
        ax.set_title(f"{side}: {branch_name} ({' -> '.join(e.name for e in branches[branch_name].elements)}) "
                     f"= {100 * q_purge[k] / q_prop[k]:.0f} % of peak propellant flow")
    axes[-1].set_xlabel("Time [s]")
    fig.suptitle(f"GN2 purge vs propellant volumetric flow (COPV at end of life: {p_copv / psi:.0f} psia, "
                 f"{t_copv:.0f} K; {regulator.name} {p_reg / psi:.0f} psia)", fontweight="bold")
    fig.tight_layout()
    fig.savefig(os.path.join(directory, "purge_flow.png"), dpi=120)
    plt.close(fig)
    return summary


if __name__ == "__main__":
    results = sys.argv[1] if len(sys.argv) > 1 else D.PROP_FILENAME
    for key, value in purge_check(results).items():
        print(f"{key:52s} {value}")
    print("wrote", os.path.join(D.PROP_RESULTS_DIRECTORY, "purge_flow.png"))
