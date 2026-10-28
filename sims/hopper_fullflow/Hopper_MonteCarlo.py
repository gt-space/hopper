"""Monte Carlo of the full prop system: every run flies the same thrust curve
with inputs drawn from MC_DISPERSIONS (HOPPER_DATA.py), in parallel across
CPU cores. Run 0 is the nominal case.

Each run is the full transient (Hopper_Prop_System.py) in its own process, so
perturbed inputs reach every module. The seed makes the draw repeatable.

MC_THRUST_CONTROL picks what the engine controller holds:
    "thrust"  the targeter is re-solved with each run's engine and MR, so the
              valves deliver whatever flow that engine needs for the target
              thrust (propellant use, residuals and Pc vary)
    "flow"    every run is commanded the nominal flows (thrust varies)

Outputs in MC_DIRECTORY:
    monte_carlo_runs.csv    one row per run: the sampled inputs and the
                            run_summary metrics (Hopper_Report.py)
    monte_carlo_stats.csv   mean, std, min, 5 / 50 / 95 %, max of each metric:
                            the error bars (e.g. LOX residual = mean +- 2 std)
    monte_carlo_envelopes.png  time histories, nominal with 5-95 % band

    python Hopper_MonteCarlo.py [runs]
"""

import csv
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np

import HOPPER_DATA as D

ENVELOPE_TRACKS = ["Thrust [N]", "Chamber Pressure [psia]", "Mixture Ratio", "COPV Pressure [psia]",
                   "LOX Tank Pressure [psia]", "Fuel Tank Pressure [psia]", "LOX Liquid Mass [kg]",
                   "Fuel Liquid Mass [kg]"]


# ---- Addressing inputs ---- #
def _elements(data):
    for group in (data.LOX_PRESSURANT_LEG, data.FUEL_PRESSURANT_LEG, data.LOX_FEED, data.FUEL_FEED):
        yield from group
    for branch in data.BRANCHES:
        yield from branch.elements
    yield from data.REGULATORS


def _resolve(data, name):
    """(object, attribute) for "GLOBAL", "GLOBAL.attr" or "<element name>.attr"."""
    if "." not in name:
        return data, name
    head, attribute = name.rsplit(".", 1)
    if hasattr(data, head):
        return getattr(data, head), attribute
    for element in _elements(data):
        if element.name == head:
            return element, attribute
    raise KeyError(f"MC_DISPERSIONS: no input or element called {head!r}")


def get_input(data, name):
    owner, attribute = _resolve(data, name)
    return getattr(owner, attribute)


def set_input(data, name, value):
    owner, attribute = _resolve(data, name)
    setattr(owner, attribute, value)


def all_cdas(data):
    """Input names of every CdA that carries flow in a run: pressurant legs,
    feeds (OTV / FTV max CdA), injectors, regen, reliefs, the regulator and
    any scheduled side-path valve. Lines (CdA from geometry) and branches with
    a valve shut for the whole run are left out."""
    from Backend.Hopper_Elements import GasLine, Line, ThrottleValve

    def element_inputs(elements):
        for element in elements:
            if isinstance(element, (Line, GasLine)):
                continue
            attribute = "max_cda" if isinstance(element, ThrottleValve) else "cda_value"
            if getattr(element, attribute, None):
                yield f"{element.name}.{attribute}"

    names = []
    for group in (data.LOX_PRESSURANT_LEG, data.FUEL_PRESSURANT_LEG, data.LOX_FEED, data.FUEL_FEED):
        names += element_inputs(group)
    for branch in data.BRANCHES:
        shut = any(getattr(e, "schedule", None) and set(e.schedule[1]) == {0.0} for e in branch.elements)
        if not shut:
            names += element_inputs(branch.elements)
    names += [f"{r.name}.cda_value" for r in data.REGULATORS if r.cda_value]
    return names


def dispersion_inputs(dispersions=None, data=D):
    """{group: [input names]} with MC_ALL_CDAS expanded."""
    dispersions = dispersions or data.MC_DISPERSIONS
    return {group: (all_cdas(data) if spec[3] == data.MC_ALL_CDAS else list(spec[3]))
            for group, spec in dispersions.items()}


def _limit(name, value):
    low, high = D.MC_INPUT_LIMITS.get(name.rsplit(".", 1)[-1], (None, None))
    if low is not None:
        value = max(value, low)
    if high is not None:
        value = min(value, high)
    return value


def sample(runs=None, seed=None, dispersions=None):
    """[{input name: value}] for every run; run 0 is nominal. "each" groups
    draw every input separately, "shared" groups move all their inputs with
    one draw. Normal draws are clipped at +-3 sigma; MC_INPUT_LIMITS apply last."""
    runs = runs or D.MC_RUNS
    dispersions = dispersions or D.MC_DISPERSIONS
    inputs = dispersion_inputs(dispersions)
    rng = np.random.default_rng(D.MC_SEED if seed is None else seed)
    nominal = {name: get_input(D, name) for names in inputs.values() for name in names}

    def unit(distribution, group):
        if distribution == "normal":
            return float(np.clip(rng.standard_normal(), -3.0, 3.0))
        if distribution == "uniform":
            return float(rng.uniform(-1.0, 1.0))
        raise ValueError(f"{group}: distribution must be 'normal' or 'uniform'")

    draws = [dict(nominal)]
    for _ in range(runs - 1):
        draw = {}
        for group, (distribution, size, mode, _, how) in dispersions.items():
            shared = unit(distribution, group) if how == "shared" else None
            for name in inputs[group]:
                u = shared if shared is not None else unit(distribution, group)
                delta = size * u * (nominal[name] if mode == "rel" else 1.0)
                draw[name] = _limit(name, nominal[name] + delta)
        draws.append(draw)
    return draws


# ---- One run (in a fresh process) ---- #
def run_case(index, draw, nominal_targets, directory, save_dt, control):
    """Apply the draw, run the transient, return (index, draw, summary, envelope)."""
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    import HOPPER_DATA as data
    for name, value in draw.items():
        set_input(data, name, value)
    # Imported after the draw so their "from HOPPER_DATA import *" sees it.
    from Backend.Hopper_Export import LoadTracks, _track
    from Backend.Hopper_Prop_System import HopperPropSystem
    from Backend.Hopper_Report import run_summary
    from Backend.Hopper_Thrust_Targeter import HopperThrustTargets

    targets = HopperThrustTargets(data.PROP_THRUST_CURVE) if control == "thrust" else nominal_targets
    filename = os.path.join(directory, f"run_{index:04d}.h5")
    start = time.perf_counter()
    try:
        HopperPropSystem(targets, save_dt=save_dt, filename=filename)
        summary = run_summary(filename, targets)
        summary["Wall time [s]"] = round(time.perf_counter() - start, 1)
        t, tracks = LoadTracks(filename)
        envelope = {"time": t, **{name: _track(tracks, name) for name in ENVELOPE_TRACKS}}
    except Exception as error:   # a failed run is reported, not fatal
        summary, envelope = {"Error": repr(error)[:300]}, None
    if not data.MC_KEEP_FILES and os.path.exists(filename):
        os.remove(filename)
    return index, draw, summary, envelope


# ---- Statistics and plots ---- #
def statistics(rows, metrics):
    stats = []
    for metric in metrics:
        values = np.array([row[metric] for row in rows if isinstance(row.get(metric), (int, float))], float)
        if values.size == 0:
            continue
        stats.append({"Metric": metric, "Nominal": rows[0].get(metric, ""), "Mean": values.mean(),
                      "Std": values.std(ddof=1) if values.size > 1 else 0.0, "Min": values.min(),
                      "P5": np.percentile(values, 5), "P50": np.percentile(values, 50),
                      "P95": np.percentile(values, 95), "Max": values.max(), "Runs": values.size})
    for row in stats:
        for key in ("Mean", "Std", "Min", "P5", "P50", "P95", "Max"):
            row[key] = float(f"{row[key]:.6g}")
    return stats


def plot_envelopes(envelopes, filename):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    good = [e for e in envelopes if e is not None]
    if not good:
        return
    grid = good[0]["time"]
    fig, axes = plt.subplots(4, 2, figsize=(14, 13), sharex=True)
    for ax, name in zip(axes.flat, ENVELOPE_TRACKS):
        curves = np.array([np.interp(grid, e["time"], e[name]) for e in good])
        p5, p95 = np.percentile(curves, [5, 95], axis=0)
        ax.fill_between(grid, curves.min(axis=0), curves.max(axis=0), color="#c6dbef", label="min-max")
        ax.fill_between(grid, p5, p95, color="#6baed6", label="5-95 %")
        ax.plot(grid, curves[0], color="black", linewidth=1.2, label="nominal (run 0)")
        ax.set_title(name)
        ax.grid(alpha=0.3)
    axes[0, 0].legend(fontsize=8)
    for ax in axes[-1]:
        ax.set_xlabel("Time [s]")
    fig.suptitle(f"Monte Carlo, {len(good)} runs", fontweight="bold")
    fig.tight_layout()
    fig.savefig(filename, dpi=120)
    plt.close(fig)


def monte_carlo(runs=None, workers=None, directory=None):
    """Run the Monte Carlo; returns (rows, stats)."""
    from Backend.Hopper_Thrust_Targeter import HopperThrustTargets
    directory = directory or D.MC_DIRECTORY
    workers = workers or D.MC_WORKERS
    os.makedirs(directory, exist_ok=True)
    draws = sample(runs)
    nominal_targets = HopperThrustTargets(D.PROP_THRUST_CURVE)
    print(f"Monte Carlo: {len(draws)} runs on {workers} processes, seed {D.MC_SEED}, control '{D.MC_THRUST_CONTROL}'")

    results = [None] * len(draws)
    start = time.perf_counter()
    with ProcessPoolExecutor(max_workers=workers, max_tasks_per_child=1) as pool:
        futures = [pool.submit(run_case, k, draw, nominal_targets, directory, D.MC_SAVE_DT, D.MC_THRUST_CONTROL)
                   for k, draw in enumerate(draws)]
        for done, future in enumerate(as_completed(futures), start=1):
            index, draw, summary, envelope = future.result()
            results[index] = (draw, summary, envelope)
            status = summary.get("Error", f"LOX residual {summary.get('LOX residual [kg]')} kg")
            print(f"  [{done}/{len(draws)}] run {index}: {status}  ({time.perf_counter() - start:.0f} s)")

    rows = [{"Run": k, **draw, **summary} for k, (draw, summary, _) in enumerate(results)]
    fields = list(dict.fromkeys(key for row in rows for key in row))
    with open(os.path.join(directory, "monte_carlo_runs.csv"), "w", newline="", encoding="utf-8") as file:
        writer = csv.DictWriter(file, fieldnames=fields, restval="")
        writer.writeheader()
        writer.writerows(rows)
    metrics = [key for key in fields if key not in ("Run", "Error", "Wall time [s]")]
    stats = statistics(rows, metrics)
    with open(os.path.join(directory, "monte_carlo_stats.csv"), "w", newline="", encoding="utf-8") as file:
        writer = csv.DictWriter(file, fieldnames=list(stats[0]))
        writer.writeheader()
        writer.writerows(stats)
    plot_envelopes([envelope for _, _, envelope in results], os.path.join(directory, "monte_carlo_envelopes.png"))
    failed = sum("Error" in summary for _, summary, _ in results)
    print(f"Done in {time.perf_counter() - start:.0f} s, {failed} failed. Results in {directory}/")
    return rows, stats


if __name__ == "__main__":
    monte_carlo(int(sys.argv[1]) if len(sys.argv) > 1 else None)
