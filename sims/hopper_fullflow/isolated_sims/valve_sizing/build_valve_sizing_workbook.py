"""Build Valve_Sizing.xlsx from the CSVs written by valve_sizing.py.

Run after valve_sizing.py, from any working directory:
    python build_valve_sizing_workbook.py

Writes only outputs/valve_sizing_outputs/Valve_Sizing.xlsx. It does not touch
Outputs/Hopper_Parameters.xlsx, which Hopper.py regenerates on every run.
"""
import csv
import json
from itertools import zip_longest
from pathlib import Path

from openpyxl import Workbook
from openpyxl.chart import Reference, ScatterChart, Series
from openpyxl.styles import Alignment, Border, Font, PatternFill, Side

SCRIPT_DIR = Path(__file__).resolve().parent
CONFIG = json.loads((SCRIPT_DIR / "valve_sizing_config.json").read_text(encoding="utf-8"))
OUTPUT_DIR = SCRIPT_DIR / "outputs" / CONFIG["output_subdirectory"]

NAVY = "1F4E78"
HEADER = "305E8F"
YELLOW = "FFF2CC"
TEXT = "1F2937"
NOTE = "55616D"
EDGE = Side(style="thin", color="B8C6D1")
BOX = Border(left=EDGE, right=EDGE, top=EDGE, bottom=EDGE)


def read_csv(name):
    with open(OUTPUT_DIR / name, newline="", encoding="utf-8") as f:
        return [{k: number(v) for k, v in row.items()} for row in csv.DictReader(f)]


def number(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return value


def by_valve(rows, valve):
    return next(r for r in rows if r["valve"] == valve)


def title(ws, text, subtitle=""):
    ws.sheet_view.showGridLines = False
    ws.sheet_properties.tabColor = NAVY
    ws.merge_cells("A2:H2")
    ws["A2"] = text
    ws["A2"].font = Font(name="Arial", size=15, bold=True, color=TEXT)
    if subtitle:
        ws.merge_cells("A3:H3")
        ws["A3"] = subtitle
        ws["A3"].font = Font(name="Arial", size=10, italic=True, color=NOTE)


def section(ws, row, text, last_col):
    ws.merge_cells(start_row=row, start_column=1, end_row=row, end_column=last_col)
    cell = ws.cell(row, 1, text)
    cell.font = Font(name="Arial", size=10, bold=True, color="FFFFFF")
    cell.fill = PatternFill("solid", fgColor=NAVY)
    cell.alignment = Alignment(vertical="center")


def table(ws, row, rows, number_format=None):
    """Write rows (first row is the header) starting at column A, boxed and styled."""
    for i, values in enumerate(rows):
        for j, value in enumerate(values):
            cell = ws.cell(row + i, j + 1, value)
            cell.border = BOX
            cell.alignment = Alignment(vertical="center", wrap_text=(i == 0))
            if i == 0:
                cell.font = Font(name="Arial", size=10, bold=True, color="FFFFFF")
                cell.fill = PatternFill("solid", fgColor=HEADER)
                cell.alignment = Alignment(horizontal="center", vertical="center", wrap_text=True)
            else:
                cell.font = Font(name="Arial", size=10, color=TEXT)
                if number_format and isinstance(value, float):
                    cell.number_format = number_format
    return row + len(rows)


def key_values(ws, row, rows):
    """Label / value / units / basis rows, label in bold."""
    for i, values in enumerate(rows):
        for j, value in enumerate(values):
            cell = ws.cell(row + i, j + 1, value)
            cell.border = BOX
            cell.font = Font(name="Arial", size=10, bold=(j == 0), color=TEXT)
            cell.alignment = Alignment(vertical="center")
    return row + len(rows)


def widths(ws, **cols):
    for col, width in cols.items():
        ws.column_dimensions[col].width = width


vent = read_csv("vent_summary.csv")
drain = read_csv("ofd_drain_summary.csv")
fill = read_csv("ofd_fill_summary.csv")
vent_sweep = read_csv("vent_cda_sweep.csv")
fill_sweep = read_csv("ofd_fill_cda_sweep.csv")
vent_history = read_csv("vent_history.csv")
ovnt = by_valve(vent, "OVNT")
fvnt = by_valve(vent, "FVNT")
ofd_fill = fill[0]
rating_psia = 500.0 * CONFIG["minimum_component_pressure_rating_factor"]
cold_k = 90.0 - CONFIG["minimum_temperature_margin_k"]

wb = Workbook()

# Summary
ws = wb.active
ws.title = "Summary"
title(ws, "Hopper valve sizing",
      "Simulation output for OVNT, FVNT, and OFD. Values are trade-study outputs, not purchase approvals.")
section(ws, 5, "Requirements and screening criteria", 4)
key_values(ws, 6, [
    ["Minimum component pressure rating", rating_psia, "psia",
     f"500 psia tank MEOP x {CONFIG['minimum_component_pressure_rating_factor']:g} rating factor"],
    ["Minimum temperature rating", cold_k, "K or colder", "LOX saturation temperature with margin"],
    ["Material compatibility", "LOX / GN2 compatible", "", "Select per vendor material documentation"],
    ["Connection / line limit", "1/2 in OD tubing", "", "Existing LOX line after QD"],
    ["OVNT / FVNT vent endpoint", CONFIG["vent_end_gauge_psig"], "psig", "Target tank pressure after vent"],
    ["OFD drain endpoint", 1, "% liquid residual", "Drain screening endpoint"],
    ["OFD fill endpoint", 100 * CONFIG["fill_end_liquid_fraction"], "% liquid fill", "Trade-study endpoint"],
    ["OVNT EOL ullage", ovnt["EOL_ullage_fraction"], "fraction", "From EOL liquid residual, not 10% loaded ullage"],
    ["FVNT EOL ullage", fvnt["EOL_ullage_fraction"], "fraction", "From EOL liquid residual, not 10% loaded ullage"],
])
ws["B13"].number_format = "0.0%"
ws["B14"].number_format = "0.0%"
section(ws, 16, "Results", 8)
table(ws, 17, [
    ["Component", "Function", "CdA [mm²]", "Result", "Result units", "Basis", "Decision use", "Source"],
    ["OVNT", "LOX ullage vent", ovnt["configured_CdA_mm2"], ovnt["vent_time_to_1psig_s"], "s",
     "EOL blowdown", "Compare valve size vs vent time", "Results tab"],
    ["FVNT", "Fuel ullage vent", fvnt["configured_CdA_mm2"], fvnt["vent_time_to_1psig_s"], "s",
     "EOL blowdown", "Compare valve size vs vent time", "Results tab"],
    ["OFD", "LOX fill", ofd_fill["example_OFD_CdA_mm2"], ofd_fill["liquid_transfer_time_s"], "s liquid transfer",
     f"{ofd_fill['LFS_supply_gauge_pressure_psig']:.0f} psig LFS, fixed restrictions",
     "Compare valve size vs fill time", "Results tab"],
], "0.0")
widths(ws, A=36, B=24, C=16, D=16, E=16, F=30, G=30, H=16)

# Inputs
ws = wb.create_sheet("Inputs")
title(ws, "Valve sizing inputs",
      "Edit valve_sizing_config.json; this sheet records the values used for this run.")
inputs = [
    ["Time step", CONFIG["time_step_s"], "s"],
    ["Vent target pressure", CONFIG["vent_end_gauge_psig"], "psig"],
    ["Vent target time", CONFIG["vent_target_time_s"], "s"],
    ["OFD drain target time", CONFIG["drain_target_time_s"], "s"],
    ["LFS supply pressure", CONFIG["lfs_supply_gauge_pressure_psig"], "psig"],
    ["LFS path CdA", CONFIG["ground_fill_path_cda_mm2"], "mm²"],
    ["OQD CdA", CONFIG["oqd_cda_mm2"], "mm²"],
    ["Fill tank ullage pressure", CONFIG["fill_tank_ullage_gauge_pressure_psig"], "psig"],
    ["Fill start fraction", CONFIG["fill_start_liquid_fraction"], "fraction"],
    ["Fill end fraction", CONFIG["fill_end_liquid_fraction"], "fraction"],
    ["Chilldown / boiloff delay", CONFIG["chilldown_and_boiloff_delay_s"], "s"],
    ["Component pressure rating factor", CONFIG["minimum_component_pressure_rating_factor"], "factor"],
    ["Minimum temperature margin", CONFIG["minimum_temperature_margin_k"], "K"],
    ["LOX density", CONFIG["liquid_properties"]["LOX_density_kg_m3"], "kg/m³"],
    ["RP-1 density", CONFIG["liquid_properties"]["RP1_density_kg_m3"], "kg/m³"],
]
end = table(ws, 5, [["Input", "Value", "Units", "Source"]] + [r + ["valve_sizing_config.json"] for r in inputs])
for r in range(6, end):
    ws.cell(r, 2).fill = PatternFill("solid", fgColor=YELLOW)
ws.freeze_panes = "A6"
widths(ws, A=32, B=16, C=16, D=34)

# Results
ws = wb.create_sheet("Results")
title(ws, "Valve sizing results",
      "Vent results use the FullFlow EOL tank state. OFD fill and drain results are separate liquid-flow trade studies.")
vent_cols = [("valve", "Valve"), ("fluid", "Fluid"), ("tank", "Tank"), ("configured_CdA_mm2", "CdA [mm²]"),
             ("water_equivalent_Cv", "Cv"), ("EOL_liquid_mass_kg", "EOL liquid [kg]"),
             ("EOL_liquid_volume_L", "EOL liquid [L]"), ("EOL_ullage_volume_L", "EOL ullage [L]"),
             ("EOL_ullage_fraction", "EOL ullage"), ("start_pressure_psia", "Start P [psia]"),
             ("ullage_temperature_K", "Ullage T [K]"), ("vent_time_to_1psig_s", "Vent time [s]"),
             ("CdA_required_for_target_mm2", "CdA for target [mm²]"), ("meets_60s_target", "Meets target")]
drain_cols = [("valve", "Valve"), ("mode", "Mode"), ("configured_CdA_mm2", "CdA [mm²]"),
              ("water_equivalent_Cv", "Cv"), ("ullage_gauge_pressure_psig", "Ullage P [psig]"),
              ("LOX_load_kg", "LOX load [kg]"), ("drain_time_to_1pct_residual_s", "Drain time [s]"),
              ("CdA_required_for_target_mm2", "CdA for target [mm²]"), ("meets_60s_target", "Meets target")]
fill_cols = [("valve", "Valve"), ("mode", "Mode"), ("example_OFD_CdA_mm2", "OFD CdA [mm²]"),
             ("LFS_path_CdA_mm2", "LFS CdA [mm²]"), ("OQD_CdA_mm2", "OQD CdA [mm²]"),
             ("combined_path_CdA_mm2", "Combined CdA [mm²]"), ("LFS_supply_gauge_pressure_psig", "Supply P [psig]"),
             ("tank_ullage_gauge_pressure_psig", "Tank P [psig]"), ("LOX_mass_added_kg", "LOX added [kg]"),
             ("fill_start_fraction", "Start fraction"), ("fill_end_fraction", "End fraction"),
             ("chilldown_and_boiloff_delay_s", "Delay [s]"), ("liquid_transfer_time_s", "Transfer [s]"),
             ("total_fill_operation_time_s", "Total [s]")]
row = 5
for name, rows, cols in [("Vent results", [ovnt, fvnt], vent_cols),
                         ("OFD drain results", drain, drain_cols),
                         ("OFD fill results", fill, fill_cols)]:
    section(ws, row, name, len(cols))
    row = table(ws, row + 1, [[h for _, h in cols]] + [[r[k] for k, _ in cols] for r in rows], "0.000") + 2
section(ws, row, "CdA trade-study data", 6)
table(ws, row + 1, [["Vent valve", "CdA [mm²]", "Vent time [s]"]]
      + [[r["valve"], r["CdA_mm2"], r["vent_time_s"]] for r in vent_sweep])
for i, values in enumerate([["OFD CdA [mm²]", "Liquid transfer time [s]"]]
                           + [[r["OFD_CdA_mm2"], r["liquid_transfer_time_s"]] for r in fill_sweep]):
    for j, value in enumerate(values):
        cell = ws.cell(row + 1 + i, 5 + j, value)
        cell.border = BOX
        cell.font = Font(name="Arial", size=10, bold=(i == 0), color="FFFFFF" if i == 0 else TEXT)
        if i == 0:
            cell.fill = PatternFill("solid", fgColor=HEADER)
            cell.alignment = Alignment(horizontal="center", wrap_text=True)
ws.freeze_panes = "A5"
widths(ws, **{c: 15 for c in "ABCDEFGHIJKLMN"})
ws.column_dimensions["B"].width = 25

# Pressure history: about 40 points per valve, plus the last point.
ws = wb.create_sheet("Pressure history")
title(ws, "EOL tank pressure versus time",
      "Nitrogen blowdown from the EOL tank state at the configured OVNT and FVNT CdA. Full resolution is in vent_history.csv.")
series = {}
for valve in ["OVNT", "FVNT"]:
    rows = [r for r in vent_history if r["valve"] == valve]
    step = max(1, len(rows) // 40)
    picked = rows[::step]
    if picked[-1] is not rows[-1]:
        picked.append(rows[-1])
    series[valve] = picked
table(ws, 5, [["OVNT time [s]", "OVNT pressure [psia]", "FVNT time [s]", "FVNT pressure [psia]"]]
      + [[a["time_s"] if a else None, a["tank_pressure_psia"] if a else None,
          b["time_s"] if b else None, b["tank_pressure_psia"] if b else None]
         for a, b in zip_longest(series["OVNT"], series["FVNT"])], "0.00")
chart = ScatterChart()
chart.title = "EOL Tank Pressure vs Time"
chart.x_axis.title = "Time [s]"
chart.y_axis.title = "Tank pressure [psia]"
for t in [chart.title, chart.x_axis.title, chart.y_axis.title]:
    t.overlay = False
chart.x_axis.delete = False
chart.y_axis.delete = False
chart.x_axis.number_format = "0"
chart.y_axis.number_format = "0"
chart.x_axis.majorGridlines = None
chart.legend.position = "r"
for col, valve, color in [(1, "OVNT", "1F77B4"), (3, "FVNT", "D62728")]:
    n = len(series[valve])
    s = Series(Reference(ws, min_col=col + 1, min_row=6, max_row=5 + n),
               Reference(ws, min_col=col, min_row=6, max_row=5 + n), title=valve)
    s.marker.symbol = "none"
    s.smooth = False
    s.graphicalProperties.line.solidFill = color
    s.graphicalProperties.line.width = 22000
    chart.series.append(s)
chart.width = 22
chart.height = 12
ws.add_chart(chart, "F5")
ws.freeze_panes = "A6"
widths(ws, A=16, B=20, C=16, D=20)

# Assumptions
ws = wb.create_sheet("Assumptions")
title(ws, "Model assumptions", "From valve_sizing_config.json for this run.")
end = table(ws, 5, [["Assumption", "Model treatment"]]
            + [[k.replace("_", " "), v] for k, v in CONFIG["assumptions"].items()])
for r in range(6, end):
    ws.cell(r, 2).alignment = Alignment(vertical="top", wrap_text=True)
    ws.cell(r, 1).alignment = Alignment(vertical="top")
    ws.row_dimensions[r].height = 15 * max(1, len(str(ws.cell(r, 2).value)) // 95 + 1)
widths(ws, A=30, B=105)

path = OUTPUT_DIR / "Valve_Sizing.xlsx"
wb.save(path)
print(f"Wrote {path}")
