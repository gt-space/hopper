"""Tank pressure sweep: the throttle-valve CdA each tank set pressure needs to
fly the whole thrust profile.

For every point of the thrust profile (the targeter's Pc and ox / fuel flows),
the tank must push the flow through the feed's fixed restrictions (lines,
injector, regen) into the chamber; whatever pressure is left is the throttle
valve's. The valve CdA needed is then mdot / sqrt(2 rho dP_valve), and it grows
without bound as dP_valve goes to zero. The tank is taken at the bottom of its
bang-bang band (set - lower_limit), the lowest pressure it holds in flight.

This is a steady-state calculation on the element lists in HOPPER_DATA.py
(same CdAs, tube sizes and Churchill friction as the transient); it ignores
start-up, shutdown and dispersions.

Outputs in PROP_RESULTS_DIRECTORY:
    tank_pressure_sweep.csv   per tank set pressure: max OTV / FTV CdA needed
    tank_pressure_sweep.png   (a) max CdA needed vs tank set pressure, with the
                              valve limit; (b) tank set pressure each profile
                              point needs; (c) valve CdA needed along the
                              profile at the current set pressures

    python isolated_sims/Hopper_Tank_Pressure_Sweep.py
"""

import csv
import math
import os
import sys
from pathlib import Path

import numpy as np

if not __package__:   # run as a script (python isolated_sims/Hopper_Tank_Pressure_Sweep.py): make HOPPER_DATA and the Backend
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))   # package importable

import HOPPER_DATA as D
from scipy.constants import milli, pound_force, psi

from Backend.Hopper_Elements import Line, ThrottleValve, tube_area
from Backend.Hopper_Prop_System import liquid_properties

# (label, valve, tank, feed, flow key in the targets)
SIDES = (("LOX", "OTV", D.LOX_TANK, D.LOX_FEED, "ox_mass_flow"),
         ("Fuel", "FTV", D.FUEL_TANK, D.FUEL_FEED, "fuel_mass_flow"))


def churchill(reynolds, relative_roughness):
    """Darcy friction factor, laminar to turbulent (the correlation the transient uses)."""
    a = (2.457 * math.log(1.0 / ((7.0 / reynolds) ** 0.9 + 0.27 * relative_roughness))) ** 16
    b = (37530.0 / reynolds) ** 16
    return 8.0 * ((8.0 / reynolds) ** 12 + (a + b) ** -1.5) ** (1.0 / 12.0)


def fixed_losses(feed, mass_flow, fluid):
    """Pressure drop [Pa] across everything in a feed except the throttle valve."""
    rho, mu = fluid["density"], fluid["dynamic_viscosity"]
    total = 0.0
    for element in feed:
        if isinstance(element, ThrottleValve):
            continue
        if isinstance(element, Line):
            area = tube_area(element.inner_diameter)
            velocity = mass_flow / (rho * area)
            reynolds = max(rho * velocity * element.inner_diameter / mu, 1.0)
            f = churchill(reynolds, element.roughness / element.inner_diameter)
            total += f * element.length / element.inner_diameter * 0.5 * rho * velocity**2
        else:
            total += (mass_flow / element.cda_value) ** 2 / (2.0 * rho)
    return total


def side_needs(target, tank, feed, key, set_pressure):
    """(valve CdA needed [m^2] (inf if impossible), valve dP [Pa]) at one profile
    point with the tank at the bottom of its band below set_pressure."""
    fluid = liquid_properties(tank.propellant, set_pressure, tank.liquid_temperature)
    mass_flow = target[key]
    valve_dp = set_pressure - tank.lower_limit - target["chamber_pressure"] - fixed_losses(feed, mass_flow, fluid)
    if valve_dp <= 0.0:
        return math.inf, valve_dp
    return mass_flow / math.sqrt(2.0 * fluid["density"] * valve_dp), valve_dp


def minimum_set_pressure(target, tank, feed, key, valve_cda):
    """Lowest tank set pressure [Pa] that delivers one profile point through a
    valve of valve_cda (properties at the tank's current set pressure)."""
    fluid = liquid_properties(tank.propellant, tank.set_pressure, tank.liquid_temperature)
    mass_flow = target[key]
    return (target["chamber_pressure"] + fixed_losses(feed, mass_flow, fluid)
            + (mass_flow / valve_cda) ** 2 / (2.0 * fluid["density"]) + tank.lower_limit)


def tank_pressure_sweep(targets, directory=D.PROP_RESULTS_DIRECTORY, pressures=None,
                        valve_cda=D.TANK_SWEEP_VALVE_CDA, margin=D.TANK_SWEEP_VALVE_MARGIN):
    """Write tank_pressure_sweep.csv / .png; returns {metric: value} for the run summary."""
    import matplotlib.pyplot as plt

    pressures = list(pressures or D.TANK_SWEEP_PRESSURES)
    os.makedirs(directory, exist_ok=True)
    times = np.array([t["time"] for t in targets])
    thrust = np.array([t["thrust"] for t in targets])

    # Sweep: the largest CdA any profile point needs, per set pressure.
    sweep = {valve: [] for _, valve, *_ in SIDES}
    for set_pressure in pressures:
        for _, valve, tank, feed, key in SIDES:
            sweep[valve].append(max(side_needs(t, tank, feed, key, set_pressure)[0] for t in targets))

    # Along the profile: set pressure each point needs (valve full open and at the
    # margin), and the CdA it needs at the tanks' current set pressures.
    need_full, need_margin, cda_now = {}, {}, {}
    for label, valve, tank, feed, key in SIDES:
        need_full[label] = np.array([minimum_set_pressure(t, tank, feed, key, valve_cda) for t in targets])
        need_margin[label] = np.array([minimum_set_pressure(t, tank, feed, key, margin * valve_cda) for t in targets])
        cda_now[valve] = np.array([side_needs(t, tank, feed, key, tank.set_pressure)[0] for t in targets])

    with open(os.path.join(directory, "tank_pressure_sweep.csv"), "w", newline="", encoding="utf-8") as file:
        writer = csv.writer(file)
        writer.writerow(["Tank set pressure [psia]"] + [f"{valve} max CdA needed [mm2]" for _, valve, *_ in SIDES]
                        + [f"{valve} % of {valve_cda / milli**2:g} mm2 valve" for _, valve, *_ in SIDES])
        for k, p in enumerate(pressures):
            cdas = [sweep[valve][k] for _, valve, *_ in SIDES]
            writer.writerow([round(p / psi, 1)] + [round(c / milli**2, 2) if c < math.inf else "inf" for c in cdas]
                            + [round(100 * c / valve_cda, 1) if c < math.inf else "inf" for c in cdas])

    # ---- Figure ---- #
    colours = {"LOX": "tab:blue", "Fuel": "tab:red", "OTV": "tab:blue", "FTV": "tab:red"}
    fig, axes = plt.subplots(3, 1, figsize=(11, 13))
    p_psia = np.array(pressures) / psi
    cap = 3.0 * valve_cda / milli**2

    ax = axes[0]
    for _, valve, *_ in SIDES:
        y = np.array([c / milli**2 if c < math.inf else np.nan for c in sweep[valve]])
        ax.plot(p_psia, np.minimum(y, cap), "o-", ms=3, color=colours[valve], label=f"{valve} max CdA needed")
    ax.axhline(valve_cda / milli**2, color="k", lw=1.2, label=f"valve full open ({valve_cda / milli**2:g} mm²)")
    ax.axhline(margin * valve_cda / milli**2, color="k", lw=1, ls="--",
               label=f"{100 * margin:.0f} % open ({margin * valve_cda / milli**2:.1f} mm²)")
    for label, valve, *_ in SIDES:
        ax.axvline(need_margin[label].max() / psi, color=colours[valve], ls=":", lw=1)
    ax.set_ylim(0, cap)
    ax.set_xlabel("Tank set pressure [psia]  (tank held at set − band)")
    ax.set_ylabel("Throttle valve CdA needed [mm²]")
    ax.set_title("(a) Largest throttle CdA the profile needs vs tank pressure (dotted: lowest set pressure at the margin)")
    ax.grid(alpha=0.3); ax.legend(fontsize=9)

    ax = axes[1]
    for label, valve, *_ in SIDES:
        ax.plot(times, need_full[label] / psi, color=colours[label], lw=1.2, label=f"{label}, valve full open")
        ax.plot(times, need_margin[label] / psi, color=colours[label], lw=1, ls="--", label=f"{label}, valve at {100 * margin:.0f} %")
    ax.axhline(D.LOX_TANK.set_pressure / psi, color="0.4", lw=1, ls=":", label="current set pressure")
    twin = ax.twinx(); twin.plot(times, thrust / pound_force, color="0.7", lw=1); twin.set_ylabel("Thrust [lbf]", color="0.5")
    ax.set_xlabel("Profile time [s]"); ax.set_ylabel("Tank set pressure needed [psia]")
    ax.set_title("(b) Tank set pressure each profile point needs (grey: thrust)")
    ax.grid(alpha=0.3); ax.legend(fontsize=9, loc="lower right")

    ax = axes[2]
    for _, valve, tank, *_ in SIDES:
        ax.plot(times, np.minimum(cda_now[valve] / milli**2, cap), color=colours[valve], lw=1.2,
                label=f"{valve} at {tank.set_pressure / psi:.0f} psia set")
    ax.axhline(valve_cda / milli**2, color="k", lw=1.2); ax.axhline(margin * valve_cda / milli**2, color="k", lw=1, ls="--")
    ax.set_xlabel("Profile time [s]"); ax.set_ylabel("Throttle valve CdA needed [mm²]")
    ax.set_title("(c) Throttle CdA needed along the profile at the current tank pressures")
    ax.grid(alpha=0.3); ax.legend(fontsize=9)

    fig.suptitle(f"Tank pressure vs throttle valve CdA ({D.PROP_THRUST_CURVE_NAME}, MR {D.MIXTURE_RATIO:g})",
                 fontweight="bold")
    fig.tight_layout()
    fig.savefig(os.path.join(directory, "tank_pressure_sweep.png"), dpi=120)
    plt.close(fig)

    summary = {}
    for label, valve, *_ in SIDES:
        summary[f"Min {label} tank set pressure, {valve} full open [psia]"] = round(need_full[label].max() / psi, 1)
        summary[f"Min {label} tank set pressure, {valve} at {100 * margin:.0f} % [psia]"] = round(need_margin[label].max() / psi, 1)
        peak = cda_now[valve].max()
        summary[f"Max {valve} CdA needed at current set pressure [mm2]"] = round(peak / milli**2, 2) if peak < math.inf else "inf"
    return summary


if __name__ == "__main__":
    from Backend.Hopper_Thrust_Targeter import HopperThrustTargets
    for key, value in tank_pressure_sweep(HopperThrustTargets(D.PROP_THRUST_CURVE)).items():
        print(f"{key:55s} {value}")
    print("wrote", os.path.join(D.PROP_RESULTS_DIRECTORY, "tank_pressure_sweep.png"))
