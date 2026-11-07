"""Component table of the whole P&ID, generated from the element lists in hopper_config.json.

Every row comes from the model's own inputs, so the table always matches what
is simulated. The source column is the "note" written on that element in
hopper_config.json ("EDIT GUESS ..." marks a value nobody has confirmed).

    python Backend/Hopper_Hardware.py      writes component_table.csv and .md to PROP_RESULTS_DIRECTORY
"""

import csv
import os
from pathlib import Path

from scipy.constants import liter, milli, psi

import sys

if not __package__:   # run as a script (python Backend/Hopper_Hardware.py): make HOPPER_DATA and the Backend
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))   # package importable

import HOPPER_DATA as D
from Backend.Hopper_Elements import (BangBangValve, CavitatingVenturi, CheckValve, Custom, GasLine, Injector, Line,
                             Orifice, Regen, ReliefValve, ThrottleValve, Valve)

COLUMNS = ["Path", "Tag", "Type", "CdA [mm2]", "Details", "Source / note"]


def _note(item):
    return getattr(item, "note", "")


def element_type(element):
    """Short type name used in the table and on the schematic."""
    for cls, name in [(BangBangValve, "Bang-bang solenoid"), (CheckValve, "Check valve"),
                      (ReliefValve, "Relief valve"), (ThrottleValve, "Throttle valve"),
                      (Injector, "Injector"), (Regen, "Regen"), (Valve, "Valve"),
                      (CavitatingVenturi, "Cavitating venturi"),
                      (GasLine, "Gas line"), (Line, "Line"), (Orifice, "Orifice"), (Custom, "Custom")]:
        if isinstance(element, cls):
            return name
    return type(element).__name__


def element_cda(element):
    """Representative CdA [m^2]: fixed CdA, max CdA for a throttle, equivalent CdA for a gas line."""
    if isinstance(element, ThrottleValve):
        return element.max_cda
    if isinstance(element, CavitatingVenturi):
        return element.cavitating_discharge_coefficient * element.throat_area
    if isinstance(element, GasLine) or hasattr(element, "cda_value"):
        return element.cda_value
    return None


def element_details(element):
    """Geometry, set points and schedule in words."""
    parts = []
    if isinstance(element, (Line, GasLine)):
        size = (f'{element.outer_diameter / 0.0254:g} in x {element.wall_thickness / 0.0254:g} in wall, '
                if element.outer_diameter is not None else "")
        parts.append(f"{size}ID {element.inner_diameter / milli:.2f} mm, L {element.length:.3f} m")
        if isinstance(element, GasLine):
            parts.append("lumped (fully rough)" if not element.inertia else "FlowTube with inertia")
    if isinstance(element, ThrottleValve):
        parts.append("max CdA; CdA set to hit the flow target")
    if isinstance(element, ReliefValve):
        parts.append(f"crack {element.crack_pressure / psi:.0f} psig, full open {element.full_open_pressure / psi:.0f} psig")
    if isinstance(element, BangBangValve):
        parts.append("opened/closed by the tank bang-bang band")
    if isinstance(element, Regen):
        rises = sorted({rise for _, rise in element.delta_t_curve})
        parts.append(f"fuel heated +{rises[0]:g} K" if len(rises) == 1 else f"fuel heated +{rises[0]:g} to {rises[-1]:g} K with Pc")
    if getattr(element, "open_time", 0) or getattr(element, "close_time", 0):
        parts.append(f"travel: open {element.open_time:g} s, close {element.close_time:g} s")
    if isinstance(element, Valve) and not isinstance(element, (BangBangValve, CheckValve)):
        schedule = element.schedule
        if not schedule:
            parts.append("open")
        elif len(set(schedule[1])) == 1:
            parts.append("closed all run" if schedule[1][0] == 0 else f"fixed at {schedule[1][0]:g}")
        else:
            parts.append(f"scheduled {list(zip(schedule[0], schedule[1]))}")
    if getattr(element, "length", None) and isinstance(element, Orifice):
        parts.append(f"passage length {element.length:.3f} m")
    return "; ".join(parts)


def hardware_rows():
    """One row per P&ID item, in flow order."""
    model = "lumped as one volume" if D.COPV_MODEL == "lumped" else f"into a {D.COPV_HEADER_VOLUME / liter:g} L common header"
    rows = [["Pressurant", bottle.name, "COPV", "" if D.COPV_MODEL == "lumped" else f"{bottle.outlet_cda / milli**2:.2f}",
             f"{bottle.volume / liter:g} L, {bottle.pressure / psi:.0f} psia, {bottle.temperature:.0f} K {D.PRESSURANT}, {model}",
             _note(bottle)] for bottle in D.COPVS]
    for tank in (D.LOX_TANK, D.FUEL_TANK):
        rows.append([tank.name, tank.name, "Tank", "",
                     f"{tank.volume / liter:g} L, {tank.liquid_mass:.2f} kg {tank.propellant}, set "
                     f"{tank.set_pressure / psi:.0f} psia +{tank.upper_limit / psi:g}/-{tank.lower_limit / psi:g}, "
                     f"collapse factor {tank.collapse_factor:g}", _note(tank)])

    def add(path, element):
        cda = element_cda(element)
        rows.append([path, element.name, element_type(element), "" if cda is None else f"{cda / milli**2:.2f}",
                     element_details(element), _note(element)])

    for path, elements in [("LOX pressurant leg", D.LOX_PRESSURANT_LEG), ("Fuel pressurant leg", D.FUEL_PRESSURANT_LEG),
                           ("LOX feed", D.LOX_FEED), ("Fuel feed", D.FUEL_FEED)]:
        for element in elements:
            add(path, element)
    for regulator in D.REGULATORS:
        rows.append([f"{regulator.source} side path", regulator.name, "Regulator",
                     "" if not regulator.cda_value else f"{regulator.cda_value / milli**2:.2f}",
                     f"ideal, set {regulator.set_pressure / psi:.0f} psig", _note(regulator)])
    for branch in D.BRANCHES:
        for element in branch.elements:
            add(f"{branch.name}: {branch.source} -> {branch.sink}", element)
    for tag, note in D.NOT_MODELLED.items():
        rows.append(["Not in the flow model", tag, "", "", "", note])
    return rows


def write_table(directory=D.PROP_RESULTS_DIRECTORY):
    """component_table.csv and component_table.md in directory; returns their paths."""
    os.makedirs(directory, exist_ok=True)
    rows = hardware_rows()
    csv_path = os.path.join(directory, "component_table.csv")
    with open(csv_path, "w", newline="", encoding="utf-8") as file:
        writer = csv.writer(file)
        writer.writerow(COLUMNS)
        writer.writerows(rows)
    md_path = os.path.join(directory, "component_table.md")
    with open(md_path, "w", encoding="utf-8") as file:
        file.write("| " + " | ".join(COLUMNS) + " |\n|" + "---|" * len(COLUMNS) + "\n")
        for row in rows:
            file.write("| " + " | ".join(str(cell).replace("|", "/") for cell in row) + " |\n")
    return csv_path, md_path


if __name__ == "__main__":
    for path in write_table():
        print("wrote", path)
