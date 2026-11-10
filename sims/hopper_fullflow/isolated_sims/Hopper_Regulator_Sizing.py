"""Regulator sizing: the smallest regulator CdA (and Cv) that lets every purge reach
its requested flow for the whole purge sequence, at every regulator set pressure,
with everything else as configured.

Same blowdown model as Hopper_Purge_Sizing.py: from the last full run it takes the
COPV state at the shutdown request; the COPV then expands isentropically as the
purges draw on it. All purges open at their configured times; overlapping purges,
and the torch purge (open throughout, to ambient), share the regulator. The
regulator holds its set (plus its supply-pressure effect) until its CdA, choked from
the COPV, cannot pass the combined demand. Each set is lowered by the regulator's
full-open droop (the outlet a proportional regulator gives when fully open), so the
answer holds with the regulator wide open.

The usable sets lie between two limits: below the path limit no regulator is big
enough (the purge lines, valves and check valves cannot pass the request at that
pressure); above the dropout limit the COPV falls below the set before the purges
end. Between them a higher set needs a smaller regulator.

Requested flow, per purge, either
    the config requirement (purge.required_fraction_of_peak_flow of the propellant's
    peak volume flow, at the feed pressure at that peak), or
    --mass-flow: grams per second each purge must reach (one value for every purge,
    or one per purge in purge.sequence order).

Outputs in PROP_RESULTS_DIRECTORY:
    regulator_sizing.csv            minimum regulator CdA / Cv at each set pressure
    regulator_sizing_sequence.csv   the purge sequence at the configured set (or the lowest working
                                    set) and its minimum CdA
    regulator_sizing.png            (a) minimum CdA / Cv versus set, with the usable set range,
                                    the current regulator and the relief crack; (b) the sequence

    python isolated_sims/Hopper_Regulator_Sizing.py [results.h5] [options]
      --mass-flow 25 25       g/s each purge must reach (default: the config requirement)
      --set 750               design-point set [psig] (default: the configured set)
      --range 0.1 50          CdA search range [mm2]
      --points 24             set pressures in the sweep
"""

import argparse
import csv
import math
import os
import sys
from pathlib import Path

import numpy as np

if not __package__:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import HOPPER_DATA as D
from scipy.constants import milli, psi

from isolated_sims.Hopper_Purge_Sizing import _branch, _regulator, run_state, simulate

IN2 = 0.0254 ** 2
CV_PER_IN2 = 38.0   # Cv = 38 * CdA [in^2], the conversion used for every valve in the config


def cv(cda):
    """Flow coefficient Cv for a CdA [m^2]."""
    return CV_PER_IN2 * cda / IN2


def _setup(state, set_psig=None, mass_flow_g_s=None):
    names = [name for name, _, _ in D.PURGE_SEQUENCE]
    regulator = _regulator(_branch(names[0]))
    if regulator is None:
        raise SystemExit("The purges do not draw from a regulator; nothing to size.")
    set_gauge = regulator.set_pressure if set_psig is None else set_psig * psi
    required = None
    if mass_flow_g_s:
        values = list(mass_flow_g_s)
        if len(values) == 1:
            values *= len(names)
        if len(values) != len(names):
            raise SystemExit(f"--mass-flow needs 1 or {len(names)} values ({', '.join(names)}).")
        required = {name: value / 1000.0 for name, value in zip(names, values)}
    return regulator, set_gauge, required


def run_with(state, effective_set, cda, required, record=False):
    """One purge sequence with the regulator at this CdA, other elements as configured."""
    return simulate(state, effective_set, None, record=record, regulator_cda=cda, required_mass_flow=required)


def minimum_regulator_cda(state, effective_set, required, cda_range, relative_tolerance=0.002):
    """{"cda", "result", "status"}: the smallest CdA [m^2] in cda_range meeting the request.
    status "ok", "at_range_limit" (the bottom already works), "short" (even the top
    fails: path limited) or "dropout" (the COPV cannot hold the set)."""
    low, high = cda_range
    top = run_with(state, effective_set, high, required)
    if top["status"] != "ok":
        return {"cda": None, "result": top, "status": top["status"]}
    if run_with(state, effective_set, low, required)["status"] == "ok":
        return {"cda": low, "result": run_with(state, effective_set, low, required, record=True),
                "status": "at_range_limit"}
    while high / low - 1.0 > relative_tolerance:
        middle = math.sqrt(low * high)
        if run_with(state, effective_set, middle, required)["status"] == "ok":
            high = middle
        else:
            low = middle
    return {"cda": high, "result": run_with(state, effective_set, high, required, record=True), "status": "ok"}


def _bisect_set(test, low, high, tolerance=1.0 * psi):
    """Boundary set [Pa gauge] between test False (at low) and True (at high)."""
    while high - low > tolerance:
        middle = 0.5 * (low + high)
        if test(middle):
            high = middle
        else:
            low = middle
    return high


def regulator_sizing(filename=None, directory=None, set_psig=None, mass_flow_g_s=None, cda_range_mm2=(0.1, 50.0),
                     points=24):
    """Minimum regulator CdA / Cv versus regulator set pressure. Writes regulator_sizing.csv,
    regulator_sizing_sequence.csv and regulator_sizing.png; returns {metric: value}."""
    import matplotlib.pyplot as plt

    directory = directory or D.PROP_RESULTS_DIRECTORY
    os.makedirs(directory, exist_ok=True)
    state = run_state(filename)
    regulator, set_gauge, required = _setup(state, set_psig, mass_flow_g_s)
    droop = regulator.droop if regulator.cda_value else 0.0
    low_cda, high_cda = cda_range_mm2[0] * milli**2, cda_range_mm2[1] * milli**2
    current = regulator.cda_value

    def status(set_value, cda):
        return run_with(state, set_value - droop, cda, required)["status"]

    # Usable set range, with the largest regulator in the search range.
    backs = [D.PURGE_BACK_PRESSURE if D.PURGE_BACK_PRESSURE is not None else leg["feed_pressure"]
             for leg in state["legs"].values()]
    floor = max(backs) - D.AMBIENT_PRESSURE + droop + 1.0 * psi
    ceiling = state["copv_pressure"] - D.AMBIENT_PRESSURE
    dropout_limit = (_bisect_set(lambda v: status(v, high_cda) == "dropout", floor, ceiling)
                     if status(ceiling, high_cda) == "dropout" else ceiling)
    below = dropout_limit - 1.0 * psi
    path_limit = (_bisect_set(lambda v: status(v, high_cda) == "ok", floor, below)
                  if status(below, high_cda) == "ok" else None)

    rows = []
    for set_value in np.linspace(floor, ceiling, points):
        answer = minimum_regulator_cda(state, set_value - droop, required, (low_cda, high_cda))
        cda = answer["cda"]
        rows.append({"Regulator set [psig]": round(set_value / psi, 1),
                     "Minimum regulator CdA [mm2]": "" if cda is None else round(cda / milli**2, 4),
                     "Minimum regulator Cv": "" if cda is None else round(cv(cda), 4),
                     "Result": {"ok": "ok", "at_range_limit": f"<= {cda_range_mm2[0]:g} mm2 works",
                                "short": "path limited", "dropout": "COPV cannot hold set"}[answer["status"]],
                     "COPV after purges [psia]": round(answer["result"]["copv_end_pressure"] / psi, 1),
                     "Gas used [g]": round(1000 * answer["result"]["gas_used"], 1)})
    with open(os.path.join(directory, "regulator_sizing.csv"), "w", newline="", encoding="utf-8") as file:
        writer = csv.DictWriter(file, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    # The configured set, and the lowest set the current regulator works at.
    at_set = minimum_regulator_cda(state, set_gauge - droop, required, (low_cda, high_cda))
    current_min_set = None
    if current and path_limit is not None and status(below, current) == "ok":
        current_min_set = _bisect_set(lambda v: status(v, current) == "ok", path_limit, below)
    # Panel (b) shows the configured set, or the lowest swept set that works if it does not.
    design_set, design_answer = set_gauge, at_set
    if at_set["cda"] is None:
        usable = [r["Regulator set [psig]"] * psi for r in rows if r["Minimum regulator CdA [mm2]"] != ""]
        if usable:
            design_set = usable[0]
            design_answer = minimum_regulator_cda(state, design_set - droop, required, (low_cda, high_cda))
    design = design_answer["result"]["rows"] if design_answer["cda"] else []
    if design:
        with open(os.path.join(directory, "regulator_sizing_sequence.csv"), "w", newline="", encoding="utf-8") as file:
            writer = csv.DictWriter(file, fieldnames=list(design[0]))
            writer.writeheader()
            writer.writerows(design)

    request = ("mass flow " + ", ".join(f"{name} {1000 * value:g} g/s" for name, value in required.items())
               if required else f">= {100 * D.PURGE_REQUIRED_FRACTION:g} % of each propellant's peak volume flow")
    summary = {"Regulator sized": regulator.name,
               "Regulator purge request": request,
               "Regulator purge sequence": ", ".join(f"{n} {d:g} s from T+{o:g}" for n, o, d in D.PURGE_SEQUENCE)
                                           + (" + torch purge throughout" if D.TORCH_LEGS else ""),
               "COPV at shutdown request [psia]": round(state["copv_pressure"] / psi, 1),
               "Regulator droop allowed for [psi]": round(droop / psi, 1),
               "Lowest usable set, any regulator size [psig]": ("none: path limited at every set" if path_limit is None
                                                                else round(path_limit / psi, 1)),
               "Highest set the COPV holds [psig]": round(dropout_limit / psi, 1)}
    label = f"at the configured {set_gauge / psi:g} psig set"
    if at_set["cda"] is None:
        summary[f"Minimum regulator CdA {label}"] = ("none: path limited, raise the set" if at_set["status"] == "short"
                                                     else "none: the COPV cannot hold this set")
    else:
        summary[f"Minimum regulator CdA {label} [mm2]"] = round(at_set["cda"] / milli**2, 3)
        summary[f"Minimum regulator Cv {label}"] = round(cv(at_set["cda"]), 4)
    if current:
        summary[f"Current {regulator.name} CdA [mm2] / Cv"] = f"{current / milli**2:.3f} / {cv(current):.3f}"
        summary[f"Lowest set that works with the current {regulator.name} [psig]"] = (
            "none: too small at every set" if current_min_set is None else round(current_min_set / psi, 1))
    relief = next((e for b in D.BRANCHES if b.source == regulator.outlet_node for e in b.elements
                   if hasattr(e, "crack_pressure")), None)

    # ---- Figure ---- #
    fig, axes = plt.subplots(2, 1, figsize=(11, 10))
    ax = axes[0]
    good = [(r["Regulator set [psig]"], r["Minimum regulator CdA [mm2]"]) for r in rows
            if r["Minimum regulator CdA [mm2]"] != ""]
    if good:
        ax.plot([g[0] for g in good], [g[1] for g in good], "o-", color="tab:blue",
                label="minimum regulator CdA (left) / Cv (right)")
    if path_limit is not None:
        ax.axvspan(floor / psi, path_limit / psi, color="0.85", label="path limited: no regulator is big enough")
    else:
        ax.axvspan(floor / psi, dropout_limit / psi, color="0.85", label="path limited at every set")
    if dropout_limit < ceiling:
        ax.axvspan(dropout_limit / psi, ceiling / psi, color="mistyrose", label="COPV cannot hold the set")
    if current:
        ax.axhline(current / milli**2, color="tab:red", ls="--", lw=1.2,
                   label=f"current {regulator.name}: {current / milli**2:.2f} mm² (Cv {cv(current):.2f})")
    ax.axvline(set_gauge / psi, color="k", ls=":", lw=1, label=f"configured set {set_gauge / psi:g} psig")
    if relief is not None:
        ax.axvline(relief.crack_pressure / psi, color="tab:orange", ls="-.", lw=1,
                   label=f"{relief.name} crack {relief.crack_pressure / psi:g} psig")
    ax.set_yscale("log")
    ax.set_xlim(floor / psi, ceiling / psi)
    ax.set_xlabel("Regulator set pressure [psig]")
    ax.set_ylabel("Minimum regulator CdA [mm²]")
    second = ax.secondary_yaxis("right", functions=(lambda x: cv(x * milli**2), lambda c: c / cv(milli**2)))
    second.set_ylabel("Minimum regulator Cv")
    ax.grid(alpha=0.3, which="both")
    ax.legend(fontsize=8, loc="upper right")
    ax.set_title(f"(a) Regulator size versus set: every purge meets {request}", fontsize=10)
    ax = axes[1]
    if design:
        for name in dict.fromkeys(r["Purge"] for r in design):
            part = [r for r in design if r["Purge"] == name]
            ax.plot([r["Time [s]"] for r in part], [r["Purge / required"] for r in part], label=name)
        ax.axhline(1.0, color="k", lw=1)
        twin = ax.twinx()
        twin.plot([r["Time [s]"] for r in design], [r["Regulator outlet [psia]"] for r in design], ".",
                  color="tab:orange", ms=3, label="regulator outlet")
        twin.set_ylabel("Regulator outlet [psia]")
        twin.legend(fontsize=8, loc="lower right")
        ax.legend(fontsize=8, loc="lower left")
        which = "the configured" if design_set == set_gauge else "the lowest working swept"
        ax.set_title(f"(b) At {which} set {design_set / psi:.0f} psig and its minimum CdA "
                     f"{design_answer['cda'] / milli**2:.2f} mm² (Cv {cv(design_answer['cda']):.3f})", fontsize=10)
    else:
        ax.text(0.5, 0.5, "No regulator size works at any set",
                transform=ax.transAxes, ha="center")
        ax.set_title("(b) At the configured set", fontsize=10)
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Purge / requested")
    ax.grid(alpha=0.3)
    fig.suptitle(f"{regulator.name} sizing: set pressure and CdA", fontweight="bold")
    fig.tight_layout()
    fig.savefig(os.path.join(directory, "regulator_sizing.png"), dpi=120)
    plt.close(fig)
    return summary


if __name__ == "__main__":
    import matplotlib
    matplotlib.use("Agg")
    parser = argparse.ArgumentParser(description="Minimum regulator CdA / Cv versus set pressure for the purge request.")
    parser.add_argument("results", nargs="?", help="results file (default: the latest run)")
    parser.add_argument("--mass-flow", nargs="+", type=float, help="g/s each purge must reach")
    parser.add_argument("--set", type=float, help="design-point set [psig] (default: the configured set)")
    parser.add_argument("--range", nargs=2, type=float, default=(0.1, 50.0), help="CdA search range [mm2]")
    parser.add_argument("--points", type=int, default=24, help="set pressures in the sweep")
    arguments = parser.parse_args()
    for key, value in regulator_sizing(arguments.results, set_psig=arguments.set, mass_flow_g_s=arguments.mass_flow,
                                       cda_range_mm2=tuple(arguments.range), points=arguments.points).items():
        print(f"{key:62s} {value}")
    print("wrote", os.path.join(D.PROP_RESULTS_DIRECTORY, "regulator_sizing.png"))
