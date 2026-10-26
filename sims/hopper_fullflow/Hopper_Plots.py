"""On-request plots of a Hopper FullFlow run.

Covers every graph in the old general_fluid_network sim's output deck
("outputs needed.pptx"), plus the plots this model added, plus any single
quantity on request.

Named plots (python Hopper_Plots.py list prints this table):
  nodes              "Network Simulation: Nodes" - pressure, temperature, mass,
                     density, quality and fill level of the tanks, COPV and
                     chamber vs time (2 x 3 grid, as the old sim's)
  feed_nodes         the same grid for the nodes along both feed lines
  ox_ladder          "Pressure Ladder - Ox" - mdot, dP, CdA, Hdot, quality and
                     state of every element from the LOX tank to the chamber
  fuel_ladder        "Pressure Ladder - Fuel" - the same for the RP-1 feed
  pressurant_ladder  "Pressure Ladder - Pressurant" - the COPV-to-tank legs
  engine             "Engine Transient Performance" - Pc, thrust vs target, MR
  stiffness          "Injector Stiffness" - injector dP / Pc
  exit_pressure      "Nozzle Exit Pressure vs Throttle" - steady sweep with the
                     flow-separation limit and ambient pressure, plus this
                     run's operating points (plot exit_pressure --mr 1.25 sets
                     the sweep's mixture ratio)
  ox_stations        pressure at every LOX-feed station vs time, and OTV CdA,
                     also as a percentage of fully open
  fuel_stations      the same for the fuel feed and FTV
  pressurization     COPV and tank pressures, bang-bang states
  tank_states        gas temperatures and remaining masses
  startup            main-valve start-up, first 2 s
"all" makes every one; "deck" makes the ones in the old sim's deck.

Custom requests (any node, element, track or component attribute):
  python Hopper_Plots.py nodes "LOX Tank" "LOX Tank Ullage" COPV
  python Hopper_Plots.py ladder OTV OINJ "Regen Jacket"
  python Hopper_Plots.py track "Chamber Pressure [psia]" "Thrust [N]"
  python Hopper_Plots.py component OTV mass_flow upstream_pressure
Options for every command: --file, --out, --units E|SI, --csv (also write the
plotted numbers), --profile-time (time from the start of the thrust profile,
the old sim's clock), --tmin / --tmax, --show.

From Python:
  from Hopper_Plots import HopperRun, make_plots
  run = HopperRun("Hopper_Prop_System.h5")
  make_plots(run, ["nodes", "ox_ladder"], "Hopper_Results", csv_out=True)

Where the numbers come from
---------------------------
Hopper_Prop_System.py tracks every element's flow, pressures, CdA, state and
inlet enthalpy, and every junction's pressure. The NODES and ELEMENTS
catalogues below are generated from the element lists in HOPPER_DATA.py and
read those tracks (through FullPlot), so any element added to the input file
appears here automatically. Properties are ThermoProp's, as in the model.
Any stored component attribute can also be plotted directly ("component").

Units: "E" (default, the old sim's choice) plots pressure in psi and
temperature in deg F; "SI" uses Pa and K. Mass flow is always kg/s, CdA mm^2,
Hdot J/s.
"""

import argparse
import csv
import re
from pathlib import Path

import fullplot as fplt
import h5py
import matplotlib
import numpy as np
from scipy.constants import liter, psi

from HOPPER_DATA import *
from Hopper_Elements import GasLine, Injector, ThrottleValve, is_dynamic_line


# Simulation time at which the thrust profile starts (profile time 0).
PROFILE_START = MAIN_VALVE_OPEN_TIME + MAIN_VALVE_RAMP_TIME
# A chamber below this multiple of ambient is treated as unlit (as vehicle_sim.Engine does).
LIT_PRESSURE_RATIO = 1.5

# The flow paths, exactly as HOPPER_DATA.py lists them. Everything below is
# generated from these, so added elements appear in the plots automatically.
FEEDS = {"LOX": (LOX_TANK, LOX_FEED), "Fuel": (FUEL_TANK, FUEL_FEED)}
LEGS = {"LOX": (LOX_TANK, LOX_PRESSURANT_LEG), "Fuel": (FUEL_TANK, FUEL_PRESSURANT_LEG)}
SIDE_KEYS = {"LOX": "ox", "Fuel": "fuel"}




# ---- Reading results ---- #
def _safe_name(name):
    """Component or track name -> HDF5 group name (FullFlow's own rule)."""
    name = re.sub(r"[\\/\x00]", "_", str(name).strip())
    return re.sub(r"\s+", "_", name)


def _as_float(values):
    """HDF5 dataset -> float array. FullFlow writes '<uninitialized>' for
    attributes that were never set; those (and any other text) become NaN."""
    array = np.asarray(values)
    if array.dtype.kind in "fiub":
        return array.astype(float)
    out = np.full(array.shape, np.nan)
    for i, value in enumerate(array.ravel()):
        try:
            out.flat[i] = float(value)
        except (TypeError, ValueError):
            pass
    return out


class HopperRun:
    """One FullFlow results file, read with FullPlot: its time vector, tracks
    and component histories. Datasets are read when first asked for, then
    cached, so opening a large file is quick."""

    def __init__(self, filename=PROP_FILENAME):
        self.filename = str(filename)
        with h5py.File(self.filename, "r") as h5:   # only to list what the file holds
            self.network = next(name for name, group in h5.items()
                                if isinstance(group, h5py.Group) and "transient" in group)
            run = h5[f"{self.network}/transient/runs/base"]
            self._track_names = sorted(run["tracks"].keys())
            self._attributes = {
                key: sorted(attribute for attribute in group.keys() if attribute != "time")
                for key, group in run["components"].items() if isinstance(group, h5py.Group)
            }
        self._file = fplt.open(self.filename).at(f"{self.network}/transient/runs/base")
        self.time = _as_float(self._file.read("time"))
        self._cache = {}

    def _read(self, path):
        if path not in self._cache:
            self._cache[path] = _as_float(self._file.read(path))
        return self._cache[path]

    @property
    def tracks(self):
        """HDF5 names of every track in the file."""
        return self._track_names

    @property
    def components(self):
        """HDF5 names of every component in the file."""
        return sorted(self._attributes)

    def attributes(self, component):
        """Attribute names stored for one component."""
        key = _safe_name(component)
        if key not in self._attributes:
            raise KeyError(f"No component '{component}'. Run 'python Hopper_Plots.py list' to see them.")
        return self._attributes[key]

    def has(self, component, attribute):
        return attribute in self._attributes.get(_safe_name(component), [])

    def component(self, component, attribute):
        """History of one component attribute, as floats on the run's time vector."""
        if attribute not in self.attributes(component):
            raise KeyError(f"Component '{component}' has no attribute '{attribute}'. "
                           f"It has: {', '.join(self.attributes(component))}")
        return self._read(f"components/{_safe_name(component)}/{attribute}")

    def track(self, name):
        """History of one Network.track() quantity, by the name it was tracked with."""
        key = _safe_name(name)
        if key not in self._track_names:
            raise KeyError(f"No track '{name}'. Run 'python Hopper_Plots.py list' to see them.")
        return self._read(f"tracks/{key}")




# ---- Nodes ---- #
# Every node is a function of the run returning SI histories of pressure [Pa],
# temperature [K], mass [kg], density [kg/m^3], quality [-] (NaN when
# single-phase) and fill level [-] (liquid volume / tank volume for the tanks,
# 1 for liquid-full line junctions, 0 for gas nodes, as in the old sim). All
# values come from tracks written by Hopper_Prop_System.py; the properties are
# ThermoProp's (the model's own property library).
def _node(pressure, temperature, mass, density, quality, fill):
    return {"pressure": pressure, "temperature": temperature, "mass": mass,
            "density": density, "quality": quality, "fill_level": fill}


def _two_phase_quality(values):
    """Quality is only defined inside the two-phase dome; outside it, NaN."""
    values = np.asarray(values, dtype=float)
    return np.where((values > 0.0) & (values < 1.0), values, np.nan)


def _copv_node(run):
    p = run.track("COPV Pressure [psia]") * psi
    return _node(p, run.track("COPV Temperature [K]"), run.track("COPV Gas Mass [kg]"),
                 run.track("COPV Density [kg/m3]"), np.full_like(p, np.nan), np.zeros_like(p))


def _tank_liquid_node(label, tank):
    def node(run):
        p = run.track(f"{label} Tank Pressure [psia]") * psi
        return _node(p, np.full_like(p, tank.liquid_temperature), run.track(f"{label} Liquid Mass [kg]"),
                     run.track(f"{label} Liquid Density [kg/m3]"),
                     _two_phase_quality(run.track(f"{label} Liquid Quality [-]")),
                     run.track(f"{label} Liquid Volume [L]") * liter / tank.volume)
    return node


def _ullage_node(label):
    def node(run):
        p = run.track(f"{label} Tank Pressure [psia]") * psi
        return _node(p, run.track(f"{label} Ullage Temperature [K]"), run.track(f"{label} Ullage Gas Mass [kg]"),
                     run.track(f"{label} Ullage Density [kg/m3]"), np.full_like(p, np.nan), np.zeros_like(p))
    return node


def _liquid_junction(label, tank, node_name):
    def node(run):
        p = run.track(f"{node_name} Pressure [Pa]")
        return _node(p, np.full_like(p, tank.liquid_temperature), np.zeros_like(p),
                     run.track(f"{label} Liquid Density [kg/m3]"),
                     _two_phase_quality(run.track(f"{label} Liquid Quality [-]")), np.ones_like(p))
    return node


def _gas_junction(node_name):
    def node(run):
        p = run.track(f"{node_name} Pressure [Pa]")
        return _node(p, run.track(f"{node_name} Temperature [K]"), np.zeros_like(p),
                     run.track(f"{node_name} Density [kg/m3]"), np.full_like(p, np.nan), np.zeros_like(p))
    return node


def _chamber_node(run):
    p = run.track("Chamber Pressure [psia]") * psi
    return _node(p, run.track("Chamber Temperature [K]"), np.zeros_like(p),
                 run.component("Main Chamber Gas", "density"), np.full_like(p, np.nan), np.zeros_like(p))


def _segments(elements):
    """A pressurant leg split as the model splits it: GasLines on their own,
    neighbouring restrictions grouped (one CompressibleOrifice each)."""
    segments = []
    for element in elements:
        if is_dynamic_line(element) or not segments or is_dynamic_line(segments[-1]):
            segments.append(element if is_dynamic_line(element) else [element])
        else:
            segments[-1].append(element)
    return segments


def _leg_nodes(elements):
    """Names of the gas junctions inside a leg (between segments)."""
    return [(segment.outlet_node if is_dynamic_line(segment) else segment[-1].outlet_node)
            for segment in _segments(elements)[:-1]]


NODES = {"COPV": _copv_node}
for _label, (_tank, _) in FEEDS.items():
    NODES[_tank.name] = _tank_liquid_node(_label, _tank)
    NODES[f"{_tank.name} Ullage"] = _ullage_node(_label)
NODES["Main Chamber"] = _chamber_node
for _label, (_tank, _feed) in FEEDS.items():
    for _element in _feed[:-1]:
        NODES[_element.outlet_node] = _liquid_junction(_label, _tank, _element.outlet_node)
for _label, (_tank, _leg) in LEGS.items():
    for _name in _leg_nodes(_leg):
        NODES[_name] = _gas_junction(_name)

# Named sets of nodes. "tanks" is the old sim's "Network Simulation: Nodes" set.
NODE_GROUPS = {
    "tanks": ["COPV", LOX_TANK.name, f"{LOX_TANK.name} Ullage", FUEL_TANK.name, f"{FUEL_TANK.name} Ullage",
              "Main Chamber"],
    **{f"{SIDE_KEYS[label]}_feed": [tank.name] + [element.outlet_node for element in feed[:-1]] + ["Main Chamber"]
       for label, (tank, feed) in FEEDS.items()},
    "pressurant": ["COPV"] + [name for _, leg in LEGS.values() for name in _leg_nodes(leg)]
                  + [f"{tank.name} Ullage" for tank, _ in LEGS.values()],
    "all": list(NODES),
}


def node_data(run, name):
    """SI histories of one node (see the notes above)."""
    if name not in NODES:
        raise KeyError(f"No node '{name}'. Nodes: {', '.join(NODES)}")
    return NODES[name](run)




# ---- Elements (pressure-ladder entries) ---- #
# Each element returns SI histories of mass_flow [kg/s], dp [Pa] (upstream
# minus downstream pressure), cda [m^2], hdot [W] (mass flow x enthalpy of the
# fluid entering it, ThermoProp reference), quality [-] and state [-], from the
# per-element tracks Hopper_Prop_System.py writes for every element. State
# follows the old sim: 1 for fixed restrictions, the open/closed command for
# solenoids and scheduled valves, and for a ThrottleValve with a full-scale
# flow, target mdot / full scale (vehicle_sim.py's "mdot"-mode ThrottleValve).
# How far a ThrottleValve is actually open is its CdA, shown as a percentage of
# fully open in ox_stations / fuel_stations. Quality is NaN: no two-phase flow.
def _tracked_element(name):
    def evaluate(run):
        mass_flow = run.track(f"{name} Mass Flow [kg/s]")
        return {"mass_flow": mass_flow,
                "dp": run.track(f"{name} Upstream Pressure [Pa]") - run.track(f"{name} Downstream Pressure [Pa]"),
                "cda": run.track(f"{name} CdA [m2]"), "state": run.track(f"{name} State [-]"),
                "hdot": mass_flow * run.track(f"{name} Inlet Enthalpy [J/kg]"),
                "quality": np.full_like(mass_flow, np.nan)}
    return evaluate


def _series(cdas):
    """Effective CdA of restrictions in series, 1 / CdA^2 = sum(1 / CdA_i^2)."""
    cdas = [np.asarray(cda, dtype=float) for cda in cdas]
    with np.errstate(divide="ignore"):
        total = sum(np.where(cda > 0.0, 1.0 / np.maximum(cda, 1e-30)**2, np.inf) for cda in cdas)
    return np.where(np.isfinite(total), 1.0 / np.sqrt(total), 0.0)


def _group_member(group, member):
    """One restriction of a lumped gas group. The model solves the group as one
    orifice of their series CdA, so the group dP is split the way the old sim's
    Series did: each element takes dP_group x (CdA_series / CdA_element)^2 while
    the group is open; while it is shut, the shut element carries it all."""
    def evaluate(run):
        data = _tracked_element(member.name)(run)
        states = [run.track(f"{element.name} State [-]") for element in group]
        open_now = np.prod(states, axis=0) > 0.5
        open_cda = _series([element.cda_value for element in group])
        shut = run.track(f"{member.name} State [-]") < 0.5
        data["dp"] = np.where(open_now, data["dp"] * (open_cda / member.cda_value) ** 2,
                              np.where(shut, data["dp"], 0.0))
        return data
    return evaluate


def _feed_total(label, tank, feed):
    """The whole feed, tank to chamber (the old sim's otv_series / ftv_series
    line). dP = tank - chamber; CdA is the equivalent CdA of the whole path;
    state is 1 while the throttle valve is commanded open."""
    throttle = next(element for element in feed if isinstance(element, ThrottleValve))

    def evaluate(run):
        mass_flow = run.track(f"{feed[0].name} Mass Flow [kg/s]")
        dp = (run.track(f"{label} Tank Pressure [psia]") - run.track("Chamber Pressure [psia]")) * psi
        density = run.track(f"{label} Liquid Density [kg/m3]")
        cda = np.where(dp > 0.0, np.abs(mass_flow) / np.sqrt(2.0 * density * np.maximum(dp, 1e-9)), 0.0)
        return {"mass_flow": mass_flow, "dp": dp, "cda": cda,
                "hdot": mass_flow * run.track(f"{label} Liquid Enthalpy [J/kg]"),
                "quality": np.full_like(mass_flow, np.nan),
                "state": (run.track(f"{throttle.name} Target Mass Flow [kg/s]") > 0.0).astype(float)}
    return evaluate


def _leg_total(label, leg):
    """The whole COPV -> ullage leg: series CdA of its segments, state = the
    bang-bang command."""
    segments = _segments(leg)

    def evaluate(run):
        mass_flow = run.track(f"{label} Pressurant Mass Flow [kg/s]")
        dp = (run.track("COPV Pressure [psia]") - run.track(f"{label} Tank Pressure [psia]")) * psi
        cdas = [run.track(f"{segment.name} CdA [m2]") if is_dynamic_line(segment)
                else run.track(" + ".join(element.name for element in segment) + " Effective CdA [m2]")
                for segment in segments]
        return {"mass_flow": mass_flow, "dp": dp, "cda": _series(cdas),
                "hdot": mass_flow * run.track("COPV Enthalpy [J/kg]"),
                "quality": np.full_like(mass_flow, np.nan), "state": run.track(f"{label} Bang-Bang State")}
    return evaluate


ELEMENTS, BRANCHES = {}, {}
for _label, (_tank, _feed) in FEEDS.items():
    _total = f"{_label} Feed (tank to chamber)"
    ELEMENTS[_total] = _feed_total(_label, _tank, _feed)
    for _element in _feed:
        ELEMENTS[_element.name] = _tracked_element(_element.name)
    BRANCHES[SIDE_KEYS[_label]] = [_total] + [element.name for element in _feed]
BRANCHES["pressurant"] = []
for _label, (_tank, _leg) in LEGS.items():
    _total = f"{_label} Pressurant Leg (COPV to ullage)"
    ELEMENTS[_total] = _leg_total(_label, _leg)
    BRANCHES["pressurant"].append(_total)
    for _segment in _segments(_leg):
        for _element in ([_segment] if is_dynamic_line(_segment) else _segment):
            ELEMENTS[_element.name] = (_tracked_element(_element.name) if is_dynamic_line(_segment)
                                       else _group_member(_segment, _element))
            BRANCHES["pressurant"].append(_element.name)


def element_data(run, name):
    """SI histories of one element (see the ELEMENTS notes above)."""
    if name not in ELEMENTS:
        raise KeyError(f"No element '{name}'. Elements: {', '.join(ELEMENTS)}")
    return ELEMENTS[name](run)


def _feed_element(side, kind):
    """The first element of a type in a feed ("ox" or "fuel"), e.g. its ThrottleValve."""
    label = next(label for label, key in SIDE_KEYS.items() if key == side)
    return next(element for element in FEEDS[label][1] if isinstance(element, kind))




# ---- Plot building blocks ---- #
def _to_pressure(values, units):
    return values / psi if units == "E" else values


def _to_temperature(values, units):
    return (values - 273.15) * 9.0 / 5.0 + 32.0 if units == "E" else values


def _same(values, units):
    return values


# (data key, column label, converter, y-axis label in E units, in SI units)
NODE_PANELS = [
    ("pressure", "Pressure", _to_pressure, "Pressure [psi]", "Pressure [Pa]"),
    ("temperature", "Temperature", _to_temperature, "Temperature [F]", "Temperature [K]"),
    ("mass", "Mass", _same, "Mass [kg]", "Mass [kg]"),
    ("density", "Density", _same, "Density [kg/m³]", "Density [kg/m³]"),
    ("quality", "Quality", _same, "Quality [-]", "Quality [-]"),
    ("fill_level", "Fill level", _same, "Fill level [-]", "Fill level [-]"),
]
ELEMENT_PANELS = [
    ("mass_flow", "mdot", _same, "mdot [kg/s]", "mdot [kg/s]"),
    ("dp", "dP", _to_pressure, "dP [psi]", "dP [Pa]"),
    ("cda", "CdA", lambda values, units: values * 1e6, "CdA [mm^2]", "CdA [mm^2]"),
    ("hdot", "Hdot", _same, "Hdot [J/s]", "Hdot [J/s]"),
    ("quality", "Q", _same, "Q [0-1]", "Q [0-1]"),
    ("state", "State", _same, "State [-]", "State [-]"),
]


class Figure:
    """A finished plot: the matplotlib figure plus the numbers it shows
    (column name -> array, starting with the time column) for CSV export."""

    def __init__(self, figure, table):
        self.figure = figure
        self.table = table


def _clock(run, profile_time):
    """Time vector and its axis label."""
    if profile_time:
        return run.time - PROFILE_START, "Profile time [s]"
    return run.time, "Time [s]"


def _apply_window(axes, window):
    if window and (window[0] is not None or window[1] is not None):
        for ax in np.atleast_1d(axes).ravel():
            ax.set_xlim(left=window[0], right=window[1])


def _grid_plot(run, names, data_function, panels, title, units, profile_time, window, empty_note):
    """The old sim's 2 x 3 overlay: one panel per quantity, one line per node/element.

    Up to six lines get a legend in every panel, as the old sim's plots did;
    more get one shared legend under the grid so the data stays visible.
    """
    import matplotlib.pyplot as plt

    time, time_label = _clock(run, profile_time)
    fig, axs = plt.subplots(2, 3, figsize=(14, 7.5), sharex=True)
    table = {time_label: time}
    data = {name: data_function(run, name) for name in names}
    shared_legend = len(names) > 6
    # One colour per line, consistent across the six panels.
    palette = plt.cm.tab20.colors if len(names) > 10 else plt.rcParams["axes.prop_cycle"].by_key()["color"]
    colors = {name: palette[i % len(palette)] for i, name in enumerate(names)}

    for ax, (key, column, convert, label_e, label_si) in zip(axs.flat, panels):
        label = label_e if units == "E" else label_si
        plotted = False
        for name in names:
            values = convert(data[name][key], units)
            table[f"{name} {label}"] = values
            if np.any(np.isfinite(values)):
                ax.plot(time, values, label=name, color=colors[name])
                plotted = True
        ax.set_ylabel(label)
        ax.grid(True, alpha=0.5)
        if plotted and not shared_legend:
            ax.legend(fontsize=7)
        elif not plotted:
            ax.text(0.5, 0.5, empty_note, transform=ax.transAxes, ha="center", va="center",
                    fontsize=9, color="0.4")
    for ax in axs[-1]:
        ax.set_xlabel(time_label)
    _apply_window(axs, window)
    fig.suptitle(title, fontsize=14)
    if shared_legend:
        handles = [plt.Line2D([], [], color=colors[name]) for name in names]
        fig.legend(handles, names, loc="lower center", ncol=min(len(names), 6), fontsize=8, frameon=False)
        fig.tight_layout(rect=[0, 0.07 + 0.025 * ((len(names) - 1) // 6), 1, 0.96])
    else:
        fig.tight_layout(rect=[0, 0, 1, 0.96])
    return Figure(fig, table)


def _expand(names, groups):
    """Allow group names (e.g. "tanks", "ox") alongside individual names.
    Duplicates (a node in two groups) are kept once, in first-seen order."""
    expanded = []
    for name in names:
        for member in (groups[name] if name in groups else [name]):
            if member not in expanded:
                expanded.append(member)
    return expanded




# ---- Plots ---- #
def plot_nodes(run, names=("tanks",), title="Network Simulation: Nodes", units="E",
               profile_time=False, window=None):
    """The old sim's node grid for any nodes (or node groups, see NODE_GROUPS)."""
    return _grid_plot(run, _expand(names, NODE_GROUPS), node_data, NODE_PANELS, title, units,
                      profile_time, window, "Single-phase throughout\n(quality undefined)")


def plot_ladder(run, names=("ox",), title="Pressure Ladder", units="E", profile_time=False, window=None):
    """The old sim's connection grid (its "pressure ladder") for any elements
    (or branches, see BRANCHES)."""
    return _grid_plot(run, _expand(names, BRANCHES), element_data, ELEMENT_PANELS, title, units,
                      profile_time, window, "No two-phase flow modelled")


def plot_engine(run, units="E", profile_time=False, window=None):
    """Engine Transient Performance: Pc, thrust vs target, MR (vehicle_sim.py's figure)."""
    import matplotlib.pyplot as plt

    time, time_label = _clock(run, profile_time)
    pc = _to_pressure(run.track("Chamber Pressure [psia]") * psi, units)
    thrust, target, mr = run.track("Thrust [N]"), run.track("Target Thrust [N]"), run.track("Mixture Ratio")
    pc_label = "Chamber Pressure [psi]" if units == "E" else "Chamber Pressure [Pa]"

    fig, axs = plt.subplots(3, 1, figsize=(12, 8), sharex=True)
    axs[0].plot(time, pc, color="red")
    axs[0].set_ylabel(pc_label)
    axs[1].plot(time, thrust, color="orange", label="Simulated Thrust")
    axs[1].plot(time, target, color="black", linestyle="--", label="Target Thrust")
    axs[1].set_ylabel("Thrust [N]")
    axs[1].legend()
    axs[2].plot(time, mr, color="purple")
    axs[2].set_ylabel("Mixture Ratio")
    axs[2].set_ylim(1.0, 2.5)
    axs[2].set_xlabel(time_label)
    for ax in axs:
        ax.grid(True)
    _apply_window(axs, window)
    fig.suptitle("Engine Transient Performance")
    fig.tight_layout()
    return Figure(fig, {time_label: time, pc_label: pc, "Thrust [N]": thrust,
                        "Target Thrust [N]": target, "Mixture Ratio": mr})


def plot_stiffness(run, units="E", profile_time=False, window=None):
    """Injector Stiffness: injector dP / Pc [%], only while the chamber is lit."""
    import matplotlib.pyplot as plt

    time, time_label = _clock(run, profile_time)
    pc = run.track("Chamber Pressure [psia]")
    lit = pc > LIT_PRESSURE_RATIO * AMBIENT_PRESSURE / psi
    ox = np.where(lit, 100.0 * run.track("LOX Injector dP [psi]") / pc, np.nan)
    fuel = np.where(lit, 100.0 * run.track("Fuel Injector dP [psi]") / pc, np.nan)

    fig, ax = plt.subplots(figsize=(10, 5))
    ax.plot(time, ox, color="g", label="OINJ Stiffness")
    ax.plot(time, fuel, color="r", label="FINJ Stiffness")
    ax.set_xlabel(time_label)
    ax.set_ylabel("Stiffness (%)")
    ax.set_title("Injector Stiffness")
    ax.legend()
    ax.grid(True)
    _apply_window(ax, window)
    fig.tight_layout()
    return Figure(fig, {time_label: time, "OINJ Stiffness [%]": ox, "FINJ Stiffness [%]": fuel})


def exit_pressure_sweep(throttle_percents=np.linspace(40.0, 100.0, 15),
                        nominal_chamber_pressure=NOMINAL_CHAMBER_PRESSURE,
                        mixture_ratio=MIXTURE_RATIO, chamber_map_filename=CHAMBER_MAP_FILENAME):
    """Steady nozzle exit pressure [Pa] and thrust [N] at each throttle setting.

    Throttle is Pc as a percentage of nominal_chamber_pressure. Each point is a
    static evaluation of the same EnginePerformance block the transient uses
    (no solve needed: Pc and MR are given).
    """
    from fullflow import Network, State, SteadyState
    from Hopper_Components import EnginePerformance

    sweep = Network("Hopper Exit Pressure Sweep")
    chamber_pressure = State(nominal_chamber_pressure)
    engine = EnginePerformance("Sweep Chamber", sweep, chamber_pressure, State(mixture_ratio),
                               chamber_map_filename=chamber_map_filename)
    exit_pressure, thrust = [], []
    for throttle in throttle_percents:
        chamber_pressure.value = nominal_chamber_pressure * throttle / 100.0
        SteadyState(sweep).static_evaluate()
        exit_pressure.append(engine["nozzle"].exit_static_pressure.value)
        thrust.append(engine["thrust"].value)
    return np.asarray(throttle_percents, dtype=float), np.array(exit_pressure), np.array(thrust)


def plot_exit_pressure(run=None, units="E", profile_time=False, window=None,
                       throttle_percents=np.linspace(40.0, 100.0, 15), mixture_ratio=MIXTURE_RATIO):
    """Nozzle Exit Pressure vs Throttle (the old sim's deck figure).

    Steady sweep at mixture_ratio (default: the target MR), with the
    separation limit (NOZZLE_SEPARATION_RATIO x ambient) and ambient pressure.
    With nominal Pc 300 psia and MR 1.8 this reproduces the old deck's figure
    to within 1 %. If a run is given, its operating points (exit pressure vs Pc /
    nominal Pc while lit, at the run's own MR) are overlaid. profile_time and
    window do not apply (throttle is the x-axis).
    """
    import matplotlib.pyplot as plt

    to_unit = (lambda value: value / psi) if units == "E" else (lambda value: value)
    unit = "psia" if units == "E" else "Pa"
    throttle, exit_pressure, thrust = exit_pressure_sweep(throttle_percents, mixture_ratio=mixture_ratio)
    separation = NOZZLE_SEPARATION_RATIO * AMBIENT_PRESSURE

    fig, ax = plt.subplots(figsize=(11, 5.5))
    ax.plot(throttle, to_unit(exit_pressure), marker="o", markersize=4,
            label=f"Nozzle Exit Pressure (steady, MR {mixture_ratio:g})")
    if run is not None:
        pc = run.component("Main Chamber Nozzle", "upstream_total_pressure")
        lit = pc > LIT_PRESSURE_RATIO * AMBIENT_PRESSURE
        run_mixture_ratio = np.nanmedian(run.track("Mixture Ratio")[lit])
        ax.scatter(100.0 * pc[lit] / NOMINAL_CHAMBER_PRESSURE,
                   to_unit(run.component("Main Chamber Nozzle", "exit_static_pressure")[lit]),
                   s=6, color="tab:orange", alpha=0.5, zorder=3,
                   label=f"This run (transient, MR {run_mixture_ratio:.3g})")
    ax.axhline(to_unit(separation), color="red", linestyle="--",
               label=f"Separation limit ({NOZZLE_SEPARATION_RATIO:g} × Pamb = {to_unit(separation):.3g} {unit})")
    ax.axhline(to_unit(AMBIENT_PRESSURE), color="gray", linestyle=":", linewidth=1,
               label=f"Pamb = {to_unit(AMBIENT_PRESSURE):.3g} {unit}")
    ax.set_xlabel(f"Throttle [%]  (100 % = Pc {NOMINAL_CHAMBER_PRESSURE / psi:.0f} psia)")
    ax.set_ylabel(f"Exit Pressure [{unit}]")
    ax.set_title("Nozzle Exit Pressure vs Throttle")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="lower right")
    fig.tight_layout()
    return Figure(fig, {"Throttle [%]": throttle,
                        "Chamber Pressure [psia]": NOMINAL_CHAMBER_PRESSURE * throttle / 100.0 / psi,
                        f"Exit Pressure [{unit}]": to_unit(exit_pressure), "Thrust [N]": thrust})


def plot_stations(run, side="ox", units="E", profile_time=False, window=None):
    """Pressure at every station along one feed line vs time, and its throttle-valve
    CdA, also read as a percentage of fully open (100 % = the valve is saturated
    and can no longer hold its flow target)."""
    import matplotlib.pyplot as plt

    time, time_label = _clock(run, profile_time)
    stations = NODE_GROUPS[f"{side}_feed"]
    throttle = _feed_element(side, ThrottleValve)
    valve = throttle.name
    max_cda = throttle.max_cda * 1e6
    label = "Pressure [psia]" if units == "E" else "Pressure [Pa]"
    table = {time_label: time}

    fig, axs = plt.subplots(2, 1, figsize=(11, 7), sharex=True)
    for station in stations:
        values = _to_pressure(node_data(run, station)["pressure"], units)
        axs[0].plot(time, values, label=station)
        table[f"{station} {label}"] = values
    axs[0].set_ylabel(label)
    axs[0].legend(fontsize=8)
    cda = element_data(run, valve)["cda"] * 1e6
    axs[1].plot(time, cda, label=valve)
    axs[1].set_ylabel(f"{valve} CdA [mm^2]")
    axs[1].set_xlabel(time_label)
    opening = axs[1].secondary_yaxis("right", functions=(lambda value: 100.0 * value / max_cda,
                                                         lambda percent: percent * max_cda / 100.0))
    opening.set_ylabel(f"Opening [% of {max_cda:.1f} mm^2]")
    table[f"{valve} CdA [mm^2]"] = cda
    table[f"{valve} Opening [%]"] = 100.0 * cda / max_cda
    for ax in axs:
        ax.grid(True)
    _apply_window(axs, window)
    fig.suptitle(f"{'LOX' if side == 'ox' else 'Fuel'} Feed System: Station Pressures")
    fig.tight_layout()
    return Figure(fig, table)


def plot_pressurization(run, units="E", profile_time=False, window=None):
    """COPV pressure, tank pressures and bang-bang states."""
    import matplotlib.pyplot as plt

    time, time_label = _clock(run, profile_time)
    unit = "psia" if units == "E" else "Pa"
    series = {
        f"COPV Pressure [{unit}]": _to_pressure(run.track("COPV Pressure [psia]") * psi, units),
        f"LOX Tank Pressure [{unit}]": _to_pressure(run.track("LOX Tank Pressure [psia]") * psi, units),
        f"Fuel Tank Pressure [{unit}]": _to_pressure(run.track("Fuel Tank Pressure [psia]") * psi, units),
        "LOX Bang-Bang State": run.track("LOX Bang-Bang State"),
        "Fuel Bang-Bang State": run.track("Fuel Bang-Bang State"),
    }
    fig, axs = plt.subplots(3, 1, figsize=(11, 9), sharex=True)
    names = list(series)
    axs[0].plot(time, series[names[0]], label="COPV")
    axs[1].plot(time, series[names[1]], label="LOX tank")
    axs[1].plot(time, series[names[2]], label="Fuel tank")
    axs[2].plot(time, series["LOX Bang-Bang State"], label="LOX side")
    axs[2].plot(time, series["Fuel Bang-Bang State"], label="Fuel side")
    for ax, ylabel in zip(axs, [f"COPV pressure [{unit}]", f"Tank pressure [{unit}]", "Bang-bang state"]):
        ax.set_ylabel(ylabel)
        ax.grid(True)
        ax.legend(fontsize=8)
    axs[-1].set_xlabel(time_label)
    _apply_window(axs, window)
    fig.suptitle("Pressurization System")
    fig.tight_layout()
    return Figure(fig, {time_label: time, **series})


def plot_tank_states(run, units="E", profile_time=False, window=None):
    """Gas temperatures and remaining propellant / pressurant masses."""
    import matplotlib.pyplot as plt

    time, time_label = _clock(run, profile_time)
    unit = "F" if units == "E" else "K"
    temperatures = {
        f"COPV Temperature [{unit}]": _to_temperature(run.track("COPV Temperature [K]"), units),
        f"LOX Ullage Temperature [{unit}]": _to_temperature(run.track("LOX Ullage Temperature [K]"), units),
        f"Fuel Ullage Temperature [{unit}]": _to_temperature(run.track("Fuel Ullage Temperature [K]"), units),
    }
    masses = {name: run.track(name) for name in
              ["LOX Liquid Mass [kg]", "Fuel Liquid Mass [kg]", "COPV Gas Mass [kg]"]}
    fig, axs = plt.subplots(2, 1, figsize=(11, 7), sharex=True)
    for name, values in temperatures.items():
        axs[0].plot(time, values, label=name.rsplit(" [", 1)[0])
    for name, values in masses.items():
        axs[1].plot(time, values, label=name.rsplit(" [", 1)[0])
    axs[0].set_ylabel(f"Temperature [{unit}]")
    axs[1].set_ylabel("Mass [kg]")
    axs[1].set_xlabel(time_label)
    for ax in axs:
        ax.grid(True)
        ax.legend(fontsize=8)
    _apply_window(axs, window)
    fig.suptitle("Tank and COPV States")
    fig.tight_layout()
    return Figure(fig, {time_label: time, **temperatures, **masses})


def plot_startup(run, units="E", profile_time=False, window=None):
    """Main-valve start-up. With no window, shows the first 2 s of the run."""
    import matplotlib.pyplot as plt

    time, time_label = _clock(run, profile_time)
    if window is None:
        window = (time[0], time[0] + 2.0)
    unit = "psia" if units == "E" else "Pa"
    injectors = [_feed_element(side, Injector) for side in ("ox", "fuel")]
    pressures = {f"{injector.name} Inlet [{unit}]":
                 _to_pressure(run.track(f"{injector.name} Upstream Pressure [Pa]"), units) for injector in injectors}
    pressures[f"Main Chamber [{unit}]"] = _to_pressure(node_data(run, "Main Chamber")["pressure"], units)
    flows = {f"{injector.name} Mass Flow [kg/s]": run.track(f"{injector.name} Mass Flow [kg/s]")
             for injector in injectors}
    fig, axs = plt.subplots(2, 1, figsize=(11, 7), sharex=True)
    for name, values in pressures.items():
        axs[0].plot(time, values, label=name.rsplit(" [", 1)[0])
    for name, values in flows.items():
        axs[1].plot(time, values, label=name.rsplit(" [", 1)[0])
    axs[0].set_ylabel(f"Pressure [{unit}]")
    axs[1].set_ylabel("Mass flow [kg/s]")
    axs[1].set_xlabel(time_label)
    for ax in axs:
        ax.grid(True)
        ax.legend(fontsize=8)
    _apply_window(axs, window)
    fig.suptitle("Start-up: Main Valves")
    fig.tight_layout()
    return Figure(fig, {time_label: time, **pressures, **flows})


def plot_tracks(run, names, units="E", profile_time=False, window=None):
    """Any tracked quantities vs time, one panel each (units as tracked)."""
    import matplotlib.pyplot as plt

    time, time_label = _clock(run, profile_time)
    fig, axs = plt.subplots(len(names), 1, figsize=(11, 2.4 * len(names) + 1.2), sharex=True, squeeze=False)
    table = {time_label: time}
    for ax, name in zip(axs[:, 0], names):
        values = run.track(name)
        ax.plot(time, values)
        ax.set_ylabel(name, fontsize=8)
        ax.grid(True)
        table[name] = values
    axs[-1, 0].set_xlabel(time_label)
    _apply_window(axs, window)
    fig.tight_layout()
    return Figure(fig, table)


def plot_component(run, component, attributes, units="E", profile_time=False, window=None):
    """Any stored component attributes vs time, one panel each (SI, as stored)."""
    import matplotlib.pyplot as plt

    time, time_label = _clock(run, profile_time)
    fig, axs = plt.subplots(len(attributes), 1, figsize=(11, 2.4 * len(attributes) + 1.2), sharex=True,
                            squeeze=False)
    table = {time_label: time}
    for ax, attribute in zip(axs[:, 0], attributes):
        values = run.component(component, attribute)
        ax.plot(time, values)
        ax.set_ylabel(attribute, fontsize=8)
        ax.grid(True)
        table[f"{component} {attribute}"] = values
    axs[-1, 0].set_xlabel(time_label)
    _apply_window(axs, window)
    fig.suptitle(component)
    fig.tight_layout()
    return Figure(fig, table)




# ---- Named plots ---- #
# name: (description, function(run, units, profile_time, window) -> Figure)
PLOTS = {
    "nodes": ("Network Simulation: Nodes (tanks, ullages, COPV, chamber)",
              lambda run, **options: plot_nodes(run, ["tanks"], "Network Simulation: Nodes", **options)),
    "feed_nodes": ("Network Simulation: Nodes along both feed lines",
                   lambda run, **options: plot_nodes(run, ["ox_feed", "fuel_feed"],
                                                     "Network Simulation: Feed-Line Nodes", **options)),
    "ox_ladder": ("Pressure Ladder - Ox: LOX tank to chamber",
                  lambda run, **options: plot_ladder(run, ["ox"], "Pressure Ladder: Ox Feed System", **options)),
    "fuel_ladder": ("Pressure Ladder - Fuel: RP-1 tank to chamber",
                    lambda run, **options: plot_ladder(run, ["fuel"], "Pressure Ladder: Fuel Feed System",
                                                       **options)),
    "pressurant_ladder": ("Pressure Ladder - Pressurant: COPV to both ullages",
                          lambda run, **options: plot_ladder(run, ["pressurant"],
                                                             "Pressure Ladder: Pressurization System", **options)),
    "engine": ("Engine Transient Performance: Pc, thrust vs target, MR", plot_engine),
    "stiffness": ("Injector Stiffness", plot_stiffness),
    "exit_pressure": ("Nozzle Exit Pressure vs Throttle (steady sweep + this run)",
                      lambda run, **options: plot_exit_pressure(run, **options)),
    "ox_stations": ("LOX feed station pressures and OTV CdA / opening",
                    lambda run, **options: plot_stations(run, "ox", **options)),
    "fuel_stations": ("Fuel feed station pressures and FTV CdA / opening",
                      lambda run, **options: plot_stations(run, "fuel", **options)),
    "pressurization": ("COPV and tank pressures, bang-bang states", plot_pressurization),
    "tank_states": ("Gas temperatures and remaining masses", plot_tank_states),
    "startup": ("Main-valve start-up (first 2 s)", plot_startup),
}

# The graphs in the old sim's output deck, in deck order.
DECK_PLOTS = ["nodes", "ox_ladder", "fuel_ladder", "pressurant_ladder", "engine", "stiffness", "exit_pressure"]


def save(figure, output_dir, name, csv_out=False, show=False):
    """Write <name>.png (and <name>.csv with the plotted numbers) to output_dir."""
    import matplotlib.pyplot as plt

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    png = output_dir / f"{name}.png"
    figure.figure.savefig(png, dpi=150)
    written = [png]
    if csv_out:
        path = output_dir / f"{name}.csv"
        columns = list(figure.table)
        rows = np.column_stack([np.asarray(figure.table[column], dtype=float) for column in columns])
        with open(path, "w", newline="", encoding="utf-8") as file:
            writer = csv.writer(file)
            writer.writerow(columns)
            writer.writerows(rows.tolist())
        written.append(path)
    if not show:
        plt.close(figure.figure)
    return written


def make_plots(run, names="all", output_dir=PROP_RESULTS_DIRECTORY, units="E", csv_out=False,
               profile_time=False, window=None, show=False, plot_options=None):
    """Make named plots (see PLOTS). names: "all", "deck", one name, or a list.
    plot_options: extra arguments for single plots, by plot name, e.g.
    {"exit_pressure": {"mixture_ratio": 1.25}}."""
    import matplotlib.pyplot as plt

    plot_options = plot_options or {}
    if names == "all":
        names = list(PLOTS)
    elif names == "deck":
        names = DECK_PLOTS
    elif isinstance(names, str):
        names = [names]
    unknown = [name for name in names if name not in PLOTS]
    if unknown:
        raise KeyError(f"Unknown plot(s) {unknown}. Plots: {', '.join(PLOTS)}")

    written = []
    for name in names:
        figure = PLOTS[name][1](run, units=units, profile_time=profile_time, window=window,
                                **plot_options.get(name, {}))
        written += save(figure, output_dir, name, csv_out, show)
    if show:
        plt.show()
    return written




# ---- Command line ---- #
def _print_list(run):
    print("Named plots (python Hopper_Plots.py plot <name> ...):")
    for name, (description, _) in PLOTS.items():
        print(f"  {name:18s} {description}")
    print('  all / deck         every plot / the old sim\'s deck plots')
    print("\nNodes (python Hopper_Plots.py nodes <node or group> ...):")
    print("  " + ", ".join(NODES))
    print("  groups: " + "; ".join(f"{group} = {', '.join(members)}" for group, members in NODE_GROUPS.items()
                                    if group != "all") + "; all")
    print("\nElements (python Hopper_Plots.py ladder <element or branch> ...):")
    print("  " + ", ".join(ELEMENTS))
    print("  branches: " + "; ".join(f"{branch} = {', '.join(members)}" for branch, members in BRANCHES.items()))
    print("\nTracks (python Hopper_Plots.py track <name> ...; the name as tracked or as listed here):")
    print("  " + ", ".join(sorted(name for name in run.tracks if name != "time")))
    print("\nComponents (python Hopper_Plots.py component <component> <attribute> ...;"
          " list --attributes <component> shows its attributes):")
    print("  " + ", ".join(name.replace("_", " ") for name in run.components))


def main(argv=None):
    common = argparse.ArgumentParser(add_help=False)
    common.add_argument("--file", default=PROP_FILENAME, help="FullFlow results file (default %(default)s)")
    common.add_argument("--out", default=PROP_RESULTS_DIRECTORY, help="output folder (default %(default)s)")
    common.add_argument("--units", choices=["E", "SI"], default="E", help="E = psi / deg F (default), SI = Pa / K")
    common.add_argument("--csv", action="store_true", help="also write the plotted numbers to CSV")
    common.add_argument("--profile-time", action="store_true",
                        help="time from the start of the thrust profile (the old sim's clock)")
    common.add_argument("--tmin", type=float, help="x-axis start [s]")
    common.add_argument("--tmax", type=float, help="x-axis end [s]")
    common.add_argument("--name", help="output file name (custom plots)")
    common.add_argument("--show", action="store_true", help="also open the plot windows")

    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    commands = parser.add_subparsers(dest="command", required=True)
    listing = commands.add_parser("list", parents=[common], help="list plots, nodes, elements, tracks, components")
    listing.add_argument("--attributes", metavar="COMPONENT", help="list one component's attributes")
    plot = commands.add_parser("plot", parents=[common], help="named plots")
    plot.add_argument("names", nargs="+")
    plot.add_argument("--mr", type=float, help="mixture ratio of the exit_pressure sweep (default: MIXTURE_RATIO)")
    commands.add_parser("nodes", parents=[common], help="node grid for any nodes/groups").add_argument("names", nargs="+")
    commands.add_parser("ladder", parents=[common], help="element grid for any elements/branches").add_argument(
        "names", nargs="+")
    commands.add_parser("track", parents=[common], help="any tracked quantities").add_argument("names", nargs="+")
    component = commands.add_parser("component", parents=[common], help="any component attributes")
    component.add_argument("component")
    component.add_argument("attributes", nargs="+")
    args = parser.parse_args(argv)

    # Without --show, render off-screen (must be chosen before pyplot is imported).
    if not args.show:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    run = HopperRun(args.file)
    options = {"units": args.units, "profile_time": args.profile_time,
               "window": (args.tmin, args.tmax) if (args.tmin is not None or args.tmax is not None) else None}

    if args.command == "list":
        if args.attributes:
            print(f"{args.attributes}: {', '.join(run.attributes(args.attributes))}")
        else:
            _print_list(run)
        return

    if args.command == "plot":
        names = args.names[0] if len(args.names) == 1 and args.names[0] in ("all", "deck") else args.names
        plot_options = {"exit_pressure": {"mixture_ratio": args.mr}} if args.mr is not None else {}
        written = make_plots(run, names, args.out, args.units, args.csv, args.profile_time, options["window"],
                             args.show, plot_options)
    else:
        if args.command == "nodes":
            figure = plot_nodes(run, args.names, "Network Simulation: Nodes", **options)
        elif args.command == "ladder":
            figure = plot_ladder(run, args.names, "Pressure Ladder", **options)
        elif args.command == "track":
            figure = plot_tracks(run, args.names, **options)
        else:
            figure = plot_component(run, args.component, args.attributes, **options)
        name = args.name or {"nodes": "custom_nodes", "ladder": "custom_ladder", "track": "custom_tracks",
                             "component": f"component_{_safe_name(getattr(args, 'component', ''))}"}[args.command]
        written = save(figure, args.out, name, args.csv, args.show)
        if args.show:
            plt.show()
    for path in written:
        print(f"wrote {path}")


if __name__ == "__main__":
    try:
        main()
    except KeyError as error:
        # Unknown plot / node / element / track / component: say which, no traceback.
        raise SystemExit(f"error: {error.args[0]}")
