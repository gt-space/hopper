"""Vent blowdowns in FullFlow: one gas volume emptying through one vent valve alone.

Switched on per vent in the config's sims section: prvnt_vent, ovnt_vent and
fvnt_vent. Separate from the burn: every other valve is closed. Where the gas
starts is set by sims.vent_start_state:

    BOL  the loaded config state
           PRVNT  COPV -> ambient         all bottles (volume_l, pressure_psia, temperature_k)
           OVNT   LOX ullage -> ambient   OTANK / FTANK pressure_psia and ullage_temperature_k;
           FVNT   fuel ullage -> ambient  ullage = volume_l - liquid_mass_kg / liquid density
                                          (at the liquid temperature, as in the burn)
    EOL  the last saved point of the burn history (files.history, after shutdown and
         tail-off): COPV pressure and temperature, or tank pressure, ullage temperature
         and ullage volume

with the vent's branch from side_paths held fully open through the burn's branch
builder (build_branches: SeriesCdA + CompressibleOrifice) to fluids.ambient_pressure_psia.
Each vent is run twice, both plotted on the same pressure and temperature graphs:

    adiabatic   mass and energy balance, no heat. Pressure-enthalpy state, so the
                gas can cool into the two-phase dome and condense
    isothermal  mass balance only; gas held at its starting temperature

In the tank vents the liquid is a fixed volume below the ullage: no heat or mass
exchange with it (no boil-off), and it does not drain.

Writes <vent>_blowdown.png/.csv and the FullFlow histories
<vent>_vent_adiabatic.h5 / <vent>_vent_isothermal.h5 to the results folder.
"""

import csv
import math
import os

import h5py
import matplotlib.pyplot as plt
import numpy as np
from fullflow import Lookup, Network, Volume
from thermoprop import Fluid

import HOPPER_DATA
from HOPPER_DATA import (AMBIENT_PRESSURE, BRANCHES, COPV_MIN_TEMPERATURES, COPV_VOLUME, COPVS, FUEL_TANK,
                         LOX_TANK, PRESSURANT, VENT_START_STATE)
from Backend.Hopper_Elements import Branch, Valve
from Backend.Hopper_Prop_System import build_branches, liquid_properties
from Backend.Hopper_Solver import transient_solver

PSI = 6894.757293168361
CALLOUTS_PSIG = [15.0, 10.0, 5.0, 2.0]
STEPS = 2000
CASES = (("adiabatic", "C1"), ("isothermal", "C0"))


def vent_branch(valve):
    """The config's side path whose Valve is named valve, held fully open."""
    for branch in BRANCHES:
        names = [e.name for e in branch.elements if type(e) is Valve]
        if valve in names and branch.sink == "ambient":
            elements = [Valve(e.name, e.cda_value) if type(e) is Valve else e for e in branch.elements]
            cda = sum(e.cda_value ** -2 for e in branch.elements) ** -0.5
            return Branch(branch.name, branch.source, elements), cda
    raise ValueError(f"side_paths has no branch to ambient with a Valve named {valve!r}.")


def eol_tracks():
    """Last saved values of the burn history's tracks, and the time they are at. The path
    is read at call time: during a solving run Hopper.py points it at this run's archive
    copy (Hopper_Results/<time>/), which is only copied to Outputs/ when the run ends."""
    filename = HOPPER_DATA.PROP_FILENAME
    if not os.path.exists(filename):
        raise SystemExit(f"vent_start_state EOL needs the burn history {filename}; run Hopper.py first.")
    with h5py.File(filename, "r") as f:
        base = f[next(iter(f))]["transient/runs/base"]
        return {d.attrs["name"]: float(d[-1]) for d in base["tracks"].values()}, float(base["time"][-1])


def source_state(source):
    """Starting pressure [Pa], temperature [K], gas volume [m^3], label and minimum temperature,
    at BOL (the loaded config state) or EOL (the end of the burn history)."""
    tanks = {"LOX Ullage": ("LOX", LOX_TANK), "Fuel Ullage": ("Fuel", FUEL_TANK)}
    if source != "COPV" and source not in tanks:
        raise ValueError(f"vent source {source!r}: only COPV, LOX Ullage and Fuel Ullage are supported.")
    if VENT_START_STATE == "EOL":
        tracks, time = eol_tracks()
        if source == "COPV":
            return (tracks["COPV Pressure [psia]"] * PSI, tracks["COPV Temperature [K]"], COPV_VOLUME,
                    f"EOL ({time:.1f} s) COPV ({len(COPVS)} x {COPVS[0].volume * 1e3:g} L)",
                    min(COPV_MIN_TEMPERATURES))
        label, tank = tanks[source]
        ullage = tracks[f"{label} Ullage Volume [L]"] / 1e3
        return (tracks[f"{label} Tank Pressure [psia]"] * PSI, tracks[f"{label} Ullage Temperature [K]"], ullage,
                f"EOL ({time:.1f} s) {tank.name} ullage ({ullage * 1e3:.2f} L)", None)
    if source == "COPV":
        bottle = COPVS[0]
        return (bottle.pressure, bottle.temperature, COPV_VOLUME,
                f"BOL COPV ({len(COPVS)} x {bottle.volume * 1e3:g} L)", min(COPV_MIN_TEMPERATURES))
    _, tank = tanks[source]
    density = liquid_properties(tank.propellant, tank.set_pressure, tank.liquid_temperature)["density"]
    ullage = tank.volume - tank.liquid_mass / density
    return tank.pressure, tank.ullage_temperature, ullage, f"BOL {tank.name} ullage ({ullage * 1e3:.2f} L)", None


def run_time(pressure, temperature, volume, cda):
    """Run length [s]: about the time to vent to 2 psig, from the isothermal choked time constant."""
    gas_constant = Fluid(PRESSURANT, pressure=pressure, temperature=temperature).gas_constant
    tau = volume / (0.6847 * cda * math.sqrt(gas_constant * temperature))
    return 1.1 * tau * (math.log(pressure / AMBIENT_PRESSURE) + 3.0)


def build(valve, case):
    branch, _ = vent_branch(valve)
    pressure, temperature, volume, _, _ = source_state(branch.source)
    network = Network(f"{valve} Vent {case}")
    if case == "adiabatic":
        h0 = Fluid(PRESSURANT, pressure=pressure, temperature=temperature).enthalpy
        gas = Lookup("Vent Gas", network, Fluid, PRESSURANT, pressure=pressure, enthalpy=h0)
        store = Volume("Vent Volume", network, pressure=gas.pressure, enthalpy=gas.enthalpy,
                       temperature=gas.temperature, volume=volume, density=gas.density,
                       internal_energy=gas.internal_energy, energy_variable="enthalpy")
    else:
        gas = Lookup("Vent Gas", network, Fluid, PRESSURANT, pressure=pressure, temperature=temperature)
        store = Volume("Vent Volume", network, pressure=gas.pressure, volume=volume, density=gas.density)
    nodes = {branch.source: {"pressure": store.pressure, "fluid": gas, "gas": True, "store": "Vent Volume"}}
    outflows, _ = build_branches(network, [], [branch], nodes, {"ambient": AMBIENT_PRESSURE})
    store.mass_flow_out = outflows["Vent Volume"][0]
    network.track("Pressure [psig]", (store.pressure - AMBIENT_PRESSURE) / PSI)
    network.track("Temperature [K]", gas.temperature)
    network.track("Gas Mass [kg]", store.mass)
    network.track("Vent Mass Flow [kg/s]", outflows["Vent Volume"][0])
    return network


def read(filename):
    with h5py.File(filename, "r") as f:
        base = f[next(iter(f))]["transient/runs/base"]
        found = {"time": base["time"][:]}
        for dataset in base["tracks"].values():
            found[dataset.attrs["name"]] = dataset[:]
    return found


def saturation_onset(r):
    """Index of the first point where the gas is at its saturation temperature, or None."""
    for k, (gauge, temperature) in enumerate(zip(r["Pressure [psig]"], r["Temperature [K]"])):
        try:
            saturation = Fluid(PRESSURANT, pressure=gauge * PSI + AMBIENT_PRESSURE, quality=1.0).temperature
        except Exception:
            continue
        if temperature <= saturation + 0.01:
            return k
    return None


def crossing(time, values, level):
    k = np.argmax(values <= level)
    if values[k] > level:
        return None
    if k == 0:
        return time[0]
    return np.interp(level, [values[k], values[k - 1]], [time[k], time[k - 1]])


def simulate(valve, results_directory):
    branch, cda = vent_branch(valve)
    pressure, temperature, volume, label, t_min = source_state(branch.source)
    t_final = run_time(pressure, temperature, volume, cda)
    dt = t_final / STEPS
    print(f"{valve} vent blowdown (FullFlow): {label}, {(pressure - AMBIENT_PRESSURE) / PSI:.0f} psig, "
          f"{temperature:g} K, CdA {cda * 1e6:.4g} mm^2; {t_final:.4g} s at {dt:.3g} s steps")
    cases = {}
    for case, _ in CASES:
        filename = os.path.join(results_directory, f"{valve.lower()}_vent_{case}.h5")
        transient_solver(build(valve, case)).solve(dt=dt, t_final=t_final, filename=filename, save_dt=dt)
        cases[case] = read(filename)
        cases[case]["saturated"] = saturation_onset(cases[case])

    unit, scale = ("min", 60.0) if t_final > 120 else ("s", 1.0)
    for case, r in cases.items():
        t, p, temp = r["time"], r["Pressure [psig]"], r["Temperature [K]"]
        print(f"  {case}: ends at {t[-1]:.4g} s, {p[-1]:.1f} psig, {temp[-1]:.0f} K; coldest {temp.min():.0f} K")
        for level in CALLOUTS_PSIG:
            when = crossing(t, p, level)
            print(f"    {level:g} psig at {when / scale:.2f} {unit}" if when is not None else f"    {level:g} psig not reached")
        if t_min is not None and crossing(t, temp, t_min) is not None:
            print(f"    below the {t_min:g} K minimum at {crossing(t, temp, t_min) / scale:.2f} {unit}")
        if r["saturated"] is not None:
            k = r["saturated"]
            print(f"    gas reaches saturation (starts to condense) at {t[k] / scale:.2f} {unit}, "
                  f"{p[k]:.0f} psig, {temp[k]:.0f} K")

    with open(os.path.join(results_directory, f"{valve.lower()}_blowdown.csv"), "w", newline="",
              encoding="utf-8") as f:
        out = csv.writer(f)
        out.writerow(["case", "time_s", "pressure_psig", "temperature_k", "mass_kg", "mass_flow_g_s"])
        for case, r in cases.items():
            for row in zip(r["time"], r["Pressure [psig]"], r["Temperature [K]"], r["Gas Mass [kg]"],
                           r["Vent Mass Flow [kg/s]"] * 1e3):
                out.writerow([case] + [f"{v:.6g}" for v in row])

    fig, ax = plt.subplots(2, 1, figsize=(9, 7.5), sharex=True)
    for case, color in CASES:
        r = cases[case]
        ax[0].semilogy(r["time"] / scale, r["Pressure [psig]"], color=color, label=case)
        ax[1].plot(r["time"] / scale, r["Temperature [K]"], color=color, label=case)
        for level in CALLOUTS_PSIG:
            when = crossing(r["time"], r["Pressure [psig]"], level)
            if when is not None:
                ax[0].plot(when / scale, level, "o", color=color, ms=4)
        if r["saturated"] is not None:
            k = r["saturated"]
            for a, key in ((ax[0], "Pressure [psig]"), (ax[1], "Temperature [K]")):
                a.plot(r["time"][k] / scale, r[key][k], "x", color=color, ms=8, mew=2)
            ax[1].annotate(f"{case} gas reaches saturation ({r['Temperature [K]'][k]:.0f} K, "
                           f"{r['Pressure [psig]'][k]:.0f} psig)\nand condenses from here",
                           (r["time"][k] / scale, r["Temperature [K]"][k]), textcoords="offset points",
                           xytext=(10, 12), fontsize=8, color=color)
    lines = []
    for level in CALLOUTS_PSIG:
        times = [crossing(cases[c]["time"], cases[c]["Pressure [psig]"], level) for c, _ in CASES]
        lines.append(f"{level:>4g} psig  " + "  ".join("  --- " if w is None else f"{w / scale:6.2f}" for w in times))
    ax[0].text(0.98, 0.95, f"[{unit}]      adiab.  isoth.\n" + "\n".join(lines), transform=ax[0].transAxes,
               ha="right", va="top", family="monospace", fontsize=9,
               bbox=dict(boxstyle="round", fc="white", ec="0.6"))
    if t_min is not None:
        ax[1].axhline(t_min, color="k", ls=":", lw=1.0, label=f"minimum {t_min:g} K")
    ax[0].set_ylabel("Pressure [psig]")
    ax[1].set_ylabel("Gas temperature [K]")
    ax[1].set_xlabel(f"Time [{unit}]")
    ax[0].set_ylim(1, None)
    ax[0].legend(loc="lower left")
    ax[1].legend(loc="lower left")
    for a in ax:
        a.grid(alpha=0.3)
    ax[0].set_title(f"{valve} vent (FullFlow): {label}, CdA {cda * 1e6:.3g} mm$^2$, "
                    f"from {(pressure - AMBIENT_PRESSURE) / PSI:.0f} psig")
    fig.tight_layout()
    fig.savefig(os.path.join(results_directory, f"{valve.lower()}_blowdown.png"), dpi=150)
    plt.close(fig)


def write(results_directory, valves):
    if VENT_START_STATE == "EOL":
        _, time = eol_tracks()
        print(f"Vent EOL state from {HOPPER_DATA.PROP_FILENAME} (last saved point, {time:.2f} s)")
    for valve in valves:
        simulate(valve, results_directory)
