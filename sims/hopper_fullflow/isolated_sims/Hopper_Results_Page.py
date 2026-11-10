"""Results page for the full-vehicle run: an overview plot, a history CSV, the maxima
table and every input, as Markdown ready to paste into Notion.

Reads the latest run (Outputs/Hopper_Prop_System.h5, its run_summary.csv) and the
config it was run with (hopper_config.json). Writes to Outputs/results_page/:
    flight_overview.png          engine, valves, tanks, COPV, purges and torch in one figure
    free flight simulation.csv   the plotted quantities at every saved step
    flight_maxima.csv            max / min of every plotted quantity and when
    results_page.md              the page: plots and files to attach, maxima, run summary, inputs

    python isolated_sims/Hopper_Results_Page.py [results.h5]
"""

import csv
import json
import os
import sys
from pathlib import Path

import numpy as np

if not __package__:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import HOPPER_DATA as D
from Backend.Hopper_Plots import _safe_name

LBF = 4.4482216152605
OUT = Path(D.PROP_RESULTS_DIRECTORY) / "results_page"

# Graph -> [(quantity label, track, scale, unit)]; the overview plot draws these panels.
GRAPHS = {
    "Mass flow": [("LOX, main injector", "OINJ Mass Flow [kg/s]", 1, "kg/s"),
                  ("Fuel, main injector", "FINJ Mass Flow [kg/s]", 1, "kg/s"),
                  ("LOX, torch", "Torch LOX Mass Flow [kg/s]", 1000, "g/s"),
                  ("Fuel, torch", "Torch Fuel Mass Flow [kg/s]", 1000, "g/s")],
    "Throttle valves": [("OTV open", "OTV Opening", 100, "%"), ("FTV open", "FTV Opening", 100, "%")],
    "Thrust": [("Thrust", "Thrust [N]", 1 / LBF, "lbf"), ("Target thrust", "Target Thrust [N]", 1 / LBF, "lbf")],
    "Chamber pressure": [("Main chamber Pc", "Chamber Pressure [psia]", 1, "psia"),
                         ("Torch chamber Pc", "Torch Chamber Pressure [psia]", 1, "psia"),
                         ("Torch over main chamber", "Torch Over Main Pressure [psi]", 1, "psi")],
    "Injector dP": [("LOX (OINJ)", "LOX Injector dP [psi]", 1, "psi"), ("Fuel (FINJ)", "Fuel Injector dP [psi]", 1, "psi")],
    "Mixture ratio": [("Main chamber MR", "Mixture Ratio", 1, ""), ("Torch MR", "Torch Mixture Ratio", 1, "")],
    "Tank pressure": [("LOX tank", "LOX Tank Pressure [psia]", 1, "psia"),
                      ("Fuel tank", "Fuel Tank Pressure [psia]", 1, "psia")],
    "Propellant load": [("LOX liquid", "LOX Liquid Mass [kg]", 1, "kg"), ("Fuel liquid", "Fuel Liquid Mass [kg]", 1, "kg")],
    "Ullage temperature": [("LOX ullage", "LOX Ullage Temperature [K]", 1, "K"),
                           ("Fuel ullage", "Fuel Ullage Temperature [K]", 1, "K")],
    "COPV": [("COPV pressure", "COPV Pressure [psia]", 1, "psia"), ("COPV gas temperature", "COPV Temperature [K]", 1, "K"),
             ("COPV gas mass", "COPV Gas Mass [kg]", 1, "kg")],
    "Pressurant flow": [("LOX pressurant (OBB)", "LOX Pressurant Mass Flow [kg/s]", 1000, "g/s"),
                        ("Fuel pressurant (FBB)", "Fuel Pressurant Mass Flow [kg/s]", 1000, "g/s")],
    "Purge regulator": [("PR-REG outlet", "PR-REG Outlet Pressure [psia]", 1, "psia"),
                        ("PR-REG opening", "PR-REG Opening", 100, "%"),
                        ("PR-REG flow", "PR-REG Mass Flow [kg/s]", 1000, "g/s")],
    "Purge flow": [("LOX purge (OPV)", "LOX Purge Mass Flow [kg/s]", 1000, "g/s"),
                   ("Fuel purge (FPV)", "Fuel Purge Mass Flow [kg/s]", 1000, "g/s"),
                   ("Torch purge (TFPV)", "Torch Purge Mass Flow [kg/s]", 1000, "g/s")],
    "Regen": [("Regen outlet temperature", "Regen Jacket Outlet Temperature [K]", 1, "K"),
              ("Regen jacket dP", "Regen Jacket dP [psi]", 1, "psi")],
    "Relief flow": [("ORV (LOX tank)", "ORV Relief Mass Flow [kg/s]", 1000, "g/s"),
                    ("FRV (fuel tank)", "FRV Relief Mass Flow [kg/s]", 1000, "g/s"),
                    ("PRRV (COPV)", "PRRV Relief Mass Flow [kg/s]", 1000, "g/s"),
                    ("ORV-1 (purge reg)", "ORV-1 Relief Mass Flow [kg/s]", 1000, "g/s")],
}


def read_history(filename):
    import h5py
    with h5py.File(filename, "r") as h5:
        base = h5["Hopper_Prop_System/transient/runs/base"]
        tracks = base["tracks"]
        history = {"Time [s]": np.asarray(base["time"], dtype=float)}
        for rows in GRAPHS.values():
            for label, track, scale, unit in rows:
                key = _safe_name(track)
                if key in tracks:
                    history[f"{label} [{unit}]" if unit else label] = np.asarray(tracks[key], dtype=float) * scale
        metadata = json.loads(h5.attrs.get("hopper_metadata", "{}") or "{}")
    return history, metadata


def column(label, unit):
    return f"{label} [{unit}]" if unit else label


def fmt(value, unit):
    if unit in ("kg/s",):
        text = f"{value:.4f}"
    elif unit in ("g/s", "kg", "psi") or unit == "":
        text = f"{value:.3g}" if abs(value) < 10 else f"{value:.1f}"
    else:
        text = f"{value:.1f}"
    return f"{text} {unit}".strip()


def maxima(history):
    t = history["Time [s]"]
    rows = []
    for graph, quantities in GRAPHS.items():
        for label, _, _, unit in quantities:
            key = column(label, unit)
            if key not in history:
                continue
            values = history[key]
            i, j = int(np.argmax(values)), int(np.argmin(values))
            rows.append({"Graph": graph, "Quantity": label, "Max": fmt(values[i], unit), "At max [s]": round(t[i], 2),
                         "Min": fmt(values[j], unit), "At min [s]": round(t[j], 2)})
    return rows


def overview(history, metadata, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    t = history["Time [s]"]
    ignition, request = D.PROFILE_START_TIME, metadata.get("profile_end", t[-1])
    panels = ["Mass flow", "Throttle valves", "Thrust", "Chamber pressure", "Injector dP", "Tank pressure",
              "COPV", "Pressurant flow", "Purge regulator", "Purge flow"]
    figure, axes = plt.subplots(5, 2, figsize=(16, 19), sharex=True)
    for ax, graph in zip(axes.flat, panels):
        twin = None
        for k, (label, _, _, unit) in enumerate(GRAPHS[graph]):
            key = column(label, unit)
            if key not in history:
                continue
            other = (graph == "COPV" and unit != "psia") or (graph == "Purge regulator" and unit != "psia") \
                or (graph == "Mass flow" and unit == "g/s") or (graph == "Chamber pressure" and unit == "psi")
            if graph == "COPV" and unit == "kg":
                continue
            target = ax
            if other:
                twin = twin or ax.twinx()
                target = twin
            style = ":" if "Target" in label else ("--" if other else "-")
            target.plot(t, history[key], ls=style, lw=1.4, color=f"C{k}", label=key)
        ax.set_title(graph, fontsize=11)
        ax.grid(alpha=0.3)
        for x in (ignition, request):
            ax.axvline(x, color="0.5", lw=0.8, ls=":")
        lines, labels = ax.get_legend_handles_labels()
        if twin:
            more = twin.get_legend_handles_labels()
            lines, labels = lines + more[0], labels + more[1]
        ax.legend(lines, labels, fontsize=7, loc="best")
    for ax in axes[-1]:
        ax.set_xlabel("Time [s]   (dotted lines: main ignition, shutdown request)")
    figure.suptitle(f"Full-vehicle flight: {D.PROP_THRUST_CURVE_NAME}", fontsize=15, fontweight="bold")
    figure.tight_layout(rect=(0, 0, 1, 0.98))
    figure.savefig(path, dpi=110)
    plt.close(figure)


# ---- Inputs, from the config ---- #
def _num(value, digits=4):
    return f"{value:.{digits}g}" if isinstance(value, float) else str(value)


def describe(spec):
    """(type, value) text for one component spec."""
    kind = spec["type"]
    if kind in ("Line", "GasLine", "HeatedLine"):
        od, wall, idia = spec.get("outer_diameter_in"), spec.get("wall_thickness_in"), spec.get("inner_diameter_in")
        size = (f"{_num(od)} in OD x {_num(wall)} in wall ({od - 2 * wall:.3f} in ID)" if od and wall
                else f"{_num(idia)} in ID")
        return ("Tube" if kind == "HeatedLine" else "Gas line" if kind == "GasLine" else "Line",
                f"{size}, {_num(spec['length_in'])} in long")
    if kind == "ThrottleValve":
        return "Throttle valve", f"{_num(spec['max_cda_mm2'])} mm² full open"
    if kind == "Regen":
        curve = ", ".join(f"+{dt:g} K at {p:g} psia" for p, dt in spec["pc_psia_delta_t_k"])
        return "Regen", f"{_num(spec['cda_mm2'])} mm²; {curve}"
    if kind == "ReliefValve":
        return "Relief valve", f"{_num(spec['cda_mm2'])} mm², cracks {_num(spec['crack_pressure_psig'])} psig"
    if kind == "BangBangValve":
        return "Bang-bang solenoid", f"{_num(spec['cda_mm2'])} mm²"
    names = {"CheckValve": "Check valve", "Orifice": "Orifice", "Injector": "Injector", "Valve": "Valve"}
    return names.get(kind, kind), f"{_num(spec['cda_mm2'])} mm²"


def table(headers, rows):
    clean = lambda x: str(x).replace("|", "/").replace("\n", " ")
    lines = ["| " + " | ".join(headers) + " |", "| " + " | ".join("---" for _ in headers) + " |"]
    lines += ["| " + " | ".join(clean(x) for x in row) + " |" for row in rows]
    return "\n".join(lines)


def component_table(specs):
    rows = []
    for spec in specs:
        kind, value = describe(spec)
        rows.append([spec["name"], kind, value, spec.get("note", "")])
    return table(["Component", "Type", "Value", "Notes"], rows)


def inputs_markdown(config):
    notes = config.get("notes", {})
    sequence, engine, profile = config["sequence"], config["engine"], config["thrust_profile"]
    selected = profile["profiles"][profile["selected"]]
    md = ["# Free-Flight Sim: Inputs"]
    md += ["## Thrust profile", table(["Input", "Value", "Notes"], [
        ["Profile", profile["selected"], selected.get("note", "")],
        ["Source", selected.get("csv", "inline"), ""],
        ["Interpolation", profile["interpolation"], "step: each point's flow held until the next"],
        ["Throttle slew limit", f"{sequence['throttle_rate_limit_lbf_s']:g} lbf/s",
         notes.get("sequence.throttle_rate_limit_lbf_s", "")],
        ["Shortest throttle move", f"{sequence['step_smoothing_s']:g} s", ""],
        ["CSV resample spacing", f"{profile['csv_min_spacing_s']:g} s", ""]])]
    md += ["## Engine", table(["Input", "Value", "Notes"], [
        ["Target mixture ratio", f"{engine['mixture_ratio']:g}", notes.get("engine.mixture_ratio", "")],
        ["Throat area", f"{engine['throat_area_in2']:.4g} in²", ""],
        ["Exit area", f"{engine['exit_area_in2']:.4g} in²",
         f"Expansion ratio {engine['exit_area_in2'] / engine['throat_area_in2']:.2f}"],
        ["c* efficiency", f"{engine['cstar_efficiency']['low']:g} at {engine['cstar_efficiency']['pc_low_psia']:g} psia → "
                          f"{engine['cstar_efficiency']['high']:g} at "
                          f"{engine['cstar_efficiency']['pc_low_psia'] + engine['cstar_efficiency']['pc_span_psi']:g} psia",
         "Linear, flat outside"],
        ["Cf efficiency", f"{engine['cf_efficiency']['low']:g} at {engine['cf_efficiency']['pc_low_psia']:g} psia → "
                          f"{engine['cf_efficiency']['high']:g} at "
                          f"{engine['cf_efficiency']['pc_low_psia'] + engine['cf_efficiency']['pc_span_psi']:g} psia",
         "Linear, flat outside"],
        ["Ambient pressure", f"{config['fluids']['ambient_pressure_psia']:g} psia", ""],
        ["Propellants / pressurant", f"{config['fluids']['oxidizer']} / {config['fluids']['fuel']} / "
                                     f"{config['fluids']['pressurant']}", ""]])]
    torch_seq = config["torch"].get("sequence", {})
    ignition = sequence["main_ignition_s"]
    rows = [["OIMV (torch ox) opens", f"T+{torch_seq.get('oxidizer_open_s', float('nan')):g} s", ""],
            ["FIMV (torch fuel) opens", f"T+{torch_seq.get('fuel_open_s', float('nan')):g} s", ""],
            ["FTV (main fuel) starts opening", f"T+{sequence['main_fuel_open_s']:g} s", notes.get("sequence.main_fuel_open_s", "")],
            ["OTV (main ox) starts opening", f"T+{sequence['main_ox_open_s']:g} s", notes.get("sequence.main_ox_open_s", "")],
            ["Main ignition (profile start)", f"T+{ignition:g} s", notes.get("sequence.main_ignition_s", "")],
            ["TFPV (torch purge) opens", f"T+{ignition + torch_seq.get('purge_open_after_ignition_s', 0):g} s", ""],
            ["OIMV / FIMV close", f"T+{ignition + torch_seq.get('close_after_ignition_s', 0):g} s", ""]]
    md += ["## Start-up sequence", table(["Event", "Time", "Notes"], rows)]
    rows = [["OTV closes", f"shutdown + {sequence['ox_close_after_shutdown_s']:g} s",
             notes.get("sequence.ox_close_after_shutdown_s", "")],
            ["FTV closes", f"shutdown + {sequence['fuel_close_after_shutdown_s']:g} s",
             notes.get("sequence.fuel_close_after_shutdown_s", "")],
            ["OTV / FTV closing time", f"{sequence['main_valve_close_s']:g} s", notes.get("sequence.main_valve_close_s", "")]]
    for leg in config["purge"]["sequence"]:
        rows.append([f"{leg['branch']} open", f"shutdown + {leg['open_after_shutdown_s']:g} s for {leg['duration_s']:g} s",
                     leg.get("note", "")])
    rows.append(["TFPV closes", f"{torch_seq.get('purge_close_after_main_purges_s', 0):g} s after the last main purge", ""])
    rows.append(["Simulated after the last valve closes", f"{sequence['tail_off_s']:g} s", ""])
    md += ["## Shutdown sequence", "Shutdown request = the end of the thrust profile.",
           table(["Event", "Time", "Notes"], rows)]
    md += ["## COPVs", table(["Bottle", "Volume", "Pressure", "Temperature", "Min temperature", "Outlet CdA", "Notes"],
                             [[b["name"], f"{b['volume_l']:g} L", f"{b['pressure_psia']:g} psia", f"{b['temperature_k']:g} K",
                               f"{b['minimum_temperature_k']:g} K", f"{b['outlet_cda_mm2']:g} mm²", b.get("note", "")]
                              for b in config["bottles"]])]
    press = config["pressurization"]
    md += [table(["Input", "Value", "Notes"], [
        ["COPV model", press["copv_model"], notes.get("pressurization.copv_model", "")],
        ["Bang-bang sample period", f"{press['controller_period_s']:g} s", notes.get("pressurization.controller_period_s", "")]])]
    md += ["## Tanks", table(["Tank", "Propellant", "Volume", "Load", "Temperature", "Set pressure", "Ullage gas",
                              "Collapse factor", "Notes"],
                             [[t["name"], t["propellant"], f"{t['volume_l']:g} L", f"{t['liquid_mass_kg']:.3f} kg",
                               f"{t['liquid_temperature_k']:g} K",
                               f"{t['set_pressure_psia']:g} -{t['deadband_low_psi']:g}/+{t['deadband_high_psi']:g} psia",
                               f"{t['ullage_temperature_k']:g} K", f"{t['collapse_factor']:g}", t.get("note", "")]
                              for t in config["tanks"].values()])]
    paths = config["flow_paths"]
    md += ["## LOX pressurant leg (COPV → LOX ullage)", component_table(paths["oxidizer_pressurant"])]
    md += ["## Fuel pressurant leg (COPV → fuel ullage)", component_table(paths["fuel_pressurant"])]
    md += ["## LOX feed", component_table(paths["oxidizer_feed"])]
    md += ["## Fuel feed", component_table(paths["fuel_feed"])]
    network = config["torch"].get("network")
    if network:
        md += ["## Torch igniter",
               f"Enabled in the flight sim: {'yes' if config['torch']['enabled'] else 'no'}. "
               "Valve strokes are each valve's open / close time."]
        for title, key in (("Torch oxidizer leg", "oxidizer_leg"), ("Torch fuel leg", "fuel_leg"),
                           ("Torch purge leg", "purge_leg")):
            leg = network[key]
            md += [f"### {title} ({leg['source']} → {leg['sink']})", component_table(leg["elements"])]
        chamber, throat = network["chamber"], network["throat"]
        md += ["### Torch chamber and throat", table(["Input", "Value", "Notes"], [
            ["Torch chamber volume", f"{chamber['volume_l']:g} L", "Solved as zero volume in the flight sim"],
            ["Torch throat area", f"{throat['area_in2']:g} in²", throat.get("note", "")]])]
    regulators = config["regulators"]
    md += ["## Purge system", table(["Regulator", "Source", "Set", "CdA", "Supply effect", "Droop", "Notes"],
                                    [[r["name"], r["source"], f"{r['set_pressure_psig']:g} psig", f"{r['cda_mm2']:g} mm²",
                                      f"+{r.get('outlet_rise_psi_per_1000_psi_inlet_drop', 0):g} psi per 1000 psi drop",
                                      f"{r.get('full_open_droop_psi', 10):g} psi", r.get("note", "")]
                                     for r in regulators])]
    outlets = {f"{r['name']} Outlet" for r in regulators}
    purge_paths = [b for b in config["side_paths"] if b["source"] in outlets]
    for branch in purge_paths:
        md += [f"### {branch['name']} ({branch['source']} → {branch.get('sink', 'ambient')})",
               component_table(branch["elements"])]
    others = [b for b in config["side_paths"] if b["source"] not in outlets]
    rows = []
    for branch in others:
        for spec in branch["elements"]:
            kind, value = describe(spec)
            rows.append([branch["name"], f"{branch['source']} → {branch.get('sink', 'ambient')}", spec["name"], kind,
                         value, spec.get("note", "")])
    md += ["## Vents, reliefs and drain", table(["Path", "From → to", "Component", "Type", "Value", "Notes"], rows)]
    return "\n\n".join(md)


def run_summary_markdown(directory):
    path = Path(directory) / "run_summary.csv"
    if not path.exists():
        return ""
    with open(path, newline="", encoding="utf-8") as stream:
        rows = [row for row in csv.reader(stream)][1:]
    keep = [r for r in rows if not r[0].startswith(("LOX drain", "OVNT", "FVNT", "Min LOX tank set", "Min Fuel tank set"))]
    return "## Run summary\n\n" + table(["Metric", "Value"], keep)


def results_page(filename=None):
    OUT.mkdir(parents=True, exist_ok=True)
    filename = filename or D.PROP_FILENAME
    history, metadata = read_history(filename)
    overview(history, metadata, OUT / "flight_overview.png")
    with open(OUT / "free flight simulation.csv", "w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(list(history))
        writer.writerows(zip(*history.values()))
    rows = maxima(history)
    with open(OUT / "flight_maxima.csv", "w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    previous, table_rows = None, []
    for row in rows:   # the graph name only on its first row, as on the ramp page
        table_rows.append([row["Graph"] if row["Graph"] != previous else "", row["Quantity"], row["Max"]])
        previous = row["Graph"]
    config = json.loads(Path(D.CONFIG_PATH).read_text(encoding="utf-8-sig"))
    profile = config["thrust_profile"]["selected"]
    page = [f"# Free-Flight Profile: {profile} (full vehicle)",
            "[image: flight_overview.png]", "[image: startup_sequence.png]", "[image: shutdown_sequence.png]",
            "[file: free flight simulation.csv]",
            table(["Graph", "Quantity", "Max"], table_rows),
            run_summary_markdown(D.PROP_RESULTS_DIRECTORY),
            inputs_markdown(config)]
    (OUT / "results_page.md").write_text("\n\n".join(p for p in page if p) + "\n", encoding="utf-8")
    return OUT


if __name__ == "__main__":
    print("wrote", results_page(sys.argv[1] if len(sys.argv) > 1 else None))
