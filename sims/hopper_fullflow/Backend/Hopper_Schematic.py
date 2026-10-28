"""Schematic of the simulated P&ID, drawn from the element lists in HOPPER_DATA.py.

It shows exactly what the model contains: change an element list and the
drawing changes with it. Each box is one element (tag, type, CdA); colours
follow element type.

    python Backend/Hopper_Schematic.py     writes schematic.png to PROP_RESULTS_DIRECTORY
"""

import os

import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
from scipy.constants import liter, milli, psi

import sys
from pathlib import Path

if not __package__:   # run as a script (python Backend/Hopper_Schematic.py): make HOPPER_DATA and the Backend
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))   # package importable

import HOPPER_DATA as D
from Backend.Hopper_Hardware import element_cda, element_type

COLOURS = {
    "Bang-bang solenoid": "#9ecae1", "Valve": "#9ecae1", "Throttle valve": "#a1d99b",
    "Check valve": "#fdd0a2", "Relief valve": "#fc9272", "Orifice": "#d9d9d9", "Injector": "#bcbddc",
    "Line": "#ffffff", "Gas line": "#ffffff", "Regulator": "#fee391", "Cavitating venturi": "#c7e9c0",
    "Custom": "#f0f0f0",
}
GAS, LIQ_OX, LIQ_FUEL = "#6a51a3", "#2171b5", "#cb181d"


def _label(element):
    cda = element_cda(element)
    text = f"{element.name}\n{element_type(element)}"
    if element_type(element) in ("Line", "Gas line"):
        text = f"{element.name}\nID {element.inner_diameter / milli:.1f} mm, L {element.length:.2f} m"
    elif cda is not None:
        text += f", {cda / milli**2:.2f} mm²"
    if hasattr(element, "crack_pressure"):
        text += f"\ncracks {element.crack_pressure / psi:.0f} psig"
    schedule = getattr(element, "schedule", None)
    if schedule and len(set(schedule[1])) == 1 and schedule[1][0] == 0:
        text += "\n(closed)"
    return text


class Drawing:
    def __init__(self):
        self.fig, self.ax = plt.subplots(figsize=(22, 17))
        self.ax.set_aspect("equal")
        self.node_y = {}
        self.ax.axis("off")

    def box(self, x, y, w, h, text, colour, size=7.5, bold=False, round_=0.08):
        self.ax.add_patch(FancyBboxPatch((x - w / 2, y - h / 2), w, h, boxstyle=f"round,pad=0,rounding_size={round_}",
                                         facecolor=colour, edgecolor="black", linewidth=1.0, zorder=3))
        self.ax.text(x, y, text, ha="center", va="center", fontsize=size, zorder=4,
                     fontweight="bold" if bold else "normal")

    def pipe(self, points, colour, width=2.0, arrow=False):
        xs, ys = zip(*points)
        self.ax.plot(xs, ys, color=colour, linewidth=width, zorder=1, solid_capstyle="round")
        if arrow:
            (x0, y0), (x1, y1) = points[-2], points[-1]
            self.ax.annotate("", xy=(x1, y1), xytext=(x0, y0), zorder=2,
                             arrowprops=dict(arrowstyle="-|>", color=colour, lw=width))

    def column(self, x, y_top, elements, colour, pitch=1.15, w=3.5, h=0.8):
        """Elements top to bottom; returns the y below the last one. Records the
        position of each junction (an element's outlet node) in self.node_y."""
        y = y_top
        for element in elements:
            self.pipe([(x, y + pitch / 2), (x, y - pitch / 2)], colour)
            self.box(x, y, w, h, _label(element), COLOURS.get(element_type(element), "white"))
            self.node_y[element.outlet_node] = (x, y - pitch / 2)
            y -= pitch
        return y + pitch / 2

    def row(self, x0, y, elements, colour, direction, sink, pitch=2.35, w=2.15, h=0.85):
        """Elements left or right of (x0, y), then an arrow to the sink label."""
        x = x0 + direction * pitch * 0.75
        self.pipe([(x0, y), (x, y)], colour)
        for k, element in enumerate(elements):
            self.box(x, y, w, h, _label(element), COLOURS.get(element_type(element), "white"), size=6.8)
            if k < len(elements) - 1:
                self.pipe([(x, y), (x + direction * pitch, y)], colour)
            x += direction * pitch
        end = x - direction * pitch * 0.35
        self.pipe([(x - direction * pitch, y), (end, y)], colour, arrow=True)
        self.ax.text(end + direction * 0.12, y, sink, ha="left" if direction > 0 else "right", va="center",
                     fontsize=7, style="italic", color=colour)


def draw(filename=None):
    d = Drawing()
    ax = d.ax
    x_ox, x_fuel, y_copv = -6.0, 6.0, 15.0

    # ---- COPV and its side paths ---- #
    d.box(0, y_copv, 3.4, 1.3, f"COPV\n{D.COPV_COUNT} x {D.COPV_UNIT_VOLUME / liter:g} L {D.PRESSURANT}\n"
          f"{D.PRESSURANT_INITIAL_PRESSURE / psi:.0f} psia", "#dadaeb", size=8.5, bold=True)
    copv_branches = [b for b in D.BRANCHES if b.source == "COPV"]
    for k, branch in enumerate(copv_branches):
        y = y_copv + 0.5 - 1.0 * k
        d.row(1.7, y, branch.elements, GAS, +1, f"{branch.name} → {branch.sink}")
    for k, regulator in enumerate(D.REGULATORS):
        y_reg = y_copv + 2.4 + 2.6 * k
        d.pipe([(0, y_copv + 0.65), (0, y_reg - 0.45)], GAS)
        d.box(0, y_reg, 2.6, 0.9, f"{regulator.name}\nRegulator, set {regulator.set_pressure / psi:.0f} psig"
              + (f"\n{regulator.cda_value / milli**2:.2f} mm²" if regulator.cda_value else ""), COLOURS["Regulator"])
        fed = [b for b in D.BRANCHES if b.source == regulator.outlet_node]
        for j, branch in enumerate(fed):
            side = -1 if j % 2 == 0 else +1
            y = y_reg + 0.95 * (j // 2) - 0.45 * ((len(fed) - 1) // 2)
            x0 = side * 1.3
            d.pipe([(0, y_reg), (x0, y_reg), (x0, y)], GAS)
            d.row(x0, y, branch.elements, GAS, side, f"{branch.name} → {branch.sink}")

    # ---- Pressurant legs, tanks, feeds ---- #
    y_bus = y_copv - 1.25
    d.pipe([(0, y_copv - 0.65), (0, y_bus), (x_ox, y_bus)], GAS)
    d.pipe([(0, y_bus), (x_fuel, y_bus)], GAS)
    y_leg_top = y_bus - 0.6
    bottoms = []
    for x, leg in [(x_ox, D.LOX_PRESSURANT_LEG), (x_fuel, D.FUEL_PRESSURANT_LEG)]:
        d.pipe([(x, y_bus), (x, y_leg_top)], GAS)
        bottoms.append(d.column(x, y_leg_top - 0.1, leg, GAS))
    leg_bottoms = list(bottoms)
    y_tank_top = min(bottoms) - 0.4
    tank_h = 3.2
    y_tank = y_tank_top - tank_h / 2
    for x, tank, label, colour, feed, side, leg_bottom in [
            (x_ox, D.LOX_TANK, "LOX", LIQ_OX, D.LOX_FEED, -1, leg_bottoms[0]),
            (x_fuel, D.FUEL_TANK, "Fuel", LIQ_FUEL, D.FUEL_FEED, +1, leg_bottoms[1])]:
        d.pipe([(x, leg_bottom), (x, y_tank_top)], GAS)
        d.box(x, y_tank, 3.2, tank_h,
              f"{tank.name}\n{tank.volume / liter:g} L, {tank.liquid_mass:g} kg {tank.propellant}\n"
              f"set {tank.set_pressure / psi:.0f} psia ±{tank.upper_limit / psi:g}\n"
              f"collapse factor {tank.collapse_factor:g}", "#f7f7f7", size=8.5, bold=True, round_=0.4)
        ullage = [b for b in D.BRANCHES if b.source == f"{label} Ullage"]
        for k, branch in enumerate(ullage):
            d.row(x + side * 1.6, y_tank + 0.8 - 1.0 * k, branch.elements, GAS, side, f"{branch.name} → {branch.sink}")
        liquid = [b for b in D.BRANCHES if b.source == f"{label} Tank Liquid"]
        for k, branch in enumerate(liquid):
            d.row(x + side * 1.6, y_tank - 1.3 - 1.0 * k, branch.elements, colour, side,
                  f"{branch.name} → {branch.sink}")
        y_end = d.column(x, y_tank - tank_h / 2 - 0.75, feed, colour)
        d.pipe([(x, y_tank - tank_h / 2), (x, y_tank - tank_h / 2 - 0.2)], colour)
        bottoms.append(y_end)
        # Branches that tap a feed junction (the torch), drawn outward at that junction.
        for element in feed:
            for branch in [b for b in D.BRANCHES if b.source == element.outlet_node]:
                xj, yj = d.node_y[element.outlet_node]
                d.ax.plot([xj], [yj], "o", color=colour, markersize=5, zorder=2)
                arrow = "⇒ into" if branch.into_sink else "→"
                d.pipe([(xj, yj), (xj + side * 1.75, yj)], colour)   # out past the feed boxes
                d.row(xj + side * 1.75, yj, branch.elements, colour, side, f"{branch.name} {arrow} {branch.sink}")

    # ---- Engine ---- #
    y_engine = min(bottoms[2:]) - 1.2
    for x, colour, feed_bottom in [(x_ox, LIQ_OX, bottoms[2]), (x_fuel, LIQ_FUEL, bottoms[3])]:
        d.pipe([(x, feed_bottom), (x, y_engine), (x * 0.25, y_engine)], colour)
    d.box(0, y_engine, 3.0, 1.2, f"Main Chamber\nthroat {D.THROAT_AREA / milli**2:.0f} mm², ε {D.EXPANSION_RATIO:.2f}\n"
          f"MR {D.MIXTURE_RATIO:g}", "#fdae6b", size=8.5, bold=True)

    # ---- Purge sinks (where the purge branches land) ---- #
    purge_sinks = sorted({b.sink for b in D.BRANCHES if b.sink != "ambient" and not b.into_sink})
    ax.text(0, y_engine - 1.1, "Purge branches discharge at: " + ", ".join(purge_sinks)
            + " (back pressure only; purge gas leaves the model)", ha="center", fontsize=7.5, style="italic")

    # ---- Legend and notes ---- #
    # Legend: a strip of element colours under the engine, line colours below it.
    legend = [(n, c) for n, c in COLOURS.items() if n not in ("Gas line", "Custom")]
    for k, (name, colour) in enumerate(legend):
        d.box(-14.0 + 2.85 * k, y_engine - 2.0, 2.6, 0.42, name, colour, size=7)
    for k, (name, colour) in enumerate([("Pressurant (GN2)", GAS), ("LOX", LIQ_OX), ("Fuel", LIQ_FUEL)]):
        d.pipe([(-15.2, y_engine - 2.8 - 0.45 * k), (-14.2, y_engine - 2.8 - 0.45 * k)], colour)
        ax.text(-14.05, y_engine - 2.8 - 0.45 * k, name, va="center", fontsize=7.5)
    not_modelled = "\n".join(f"{tag}: {note}" for tag, note in D.NOT_MODELLED.items())
    ax.text(15.3, y_engine - 2.55, "On the P&ID, not in the flow model:\n" + not_modelled, ha="right", va="top",
            fontsize=7, bbox=dict(facecolor="#f7f7f7", edgecolor="#999999"))
    ax.set_title("Hopper FullFlow model schematic (generated from HOPPER_DATA.py)", fontsize=14, fontweight="bold")

    ax.autoscale_view()
    ax.margins(0.02)
    filename = filename or os.path.join(D.PROP_RESULTS_DIRECTORY, "schematic.png")
    os.makedirs(os.path.dirname(filename) or ".", exist_ok=True)
    d.fig.savefig(filename, dpi=130, bbox_inches="tight")
    plt.close(d.fig)
    return filename


if __name__ == "__main__":
    print("wrote", draw())
