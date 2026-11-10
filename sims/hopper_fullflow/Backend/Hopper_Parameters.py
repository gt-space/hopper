"""Parameter workbook: every input, component, assumption and open item in one
spreadsheet, generated from hopper_config.json so it always matches the model.

hopper_config.json stays the single source of truth: change a value there and
re-run this. The Status column is read from the note on each input ("GUESS" =
nobody has confirmed it, "EDIT" = expected to change); Owner team is a starting
suggestion by section, for whoever keeps the numbers current with the fluids,
engines and injectors teams.

    python Backend/Hopper_Parameters.py      writes Hopper_Parameters.xlsx to PROP_RESULTS_DIRECTORY
"""

import datetime
import os
from pathlib import Path

from openpyxl import Workbook
from openpyxl.styles import Alignment, Font, PatternFill
from openpyxl.utils import get_column_letter

import sys

if not __package__:   # run as a script (python Backend/Hopper_Parameters.py): make HOPPER_DATA and the Backend
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))   # package importable

import HOPPER_DATA as D
from Backend.Hopper_Config_Fields import GROUPS, UNITS
from Backend.Hopper_Hardware import COLUMNS as COMPONENT_COLUMNS, hardware_rows


OWNERS = {   # section in hopper_config.json -> team that owns those numbers
    "simulation": "Sim", "files": "Sim", "outputs": "Sim", "report": "Sim", "numerics": "Sim", "maximum_thrust": "Sim",
    "monte_carlo": "Sim", "thrust_profile": "GNC / trajectory", "fluids": "Fluids", "bottles": "Fluids", "tanks": "Fluids",
    "pressurization": "Fluids", "vent": "Fluids", "lox_drain": "Fluids", "tank_sweep": "Fluids", "torch": "Engines",
    "engine": "Engines", "chamber_map": "Engines", "sequence": "Avionics / engines",
}
SECTIONS_ON_OTHER_SHEETS = {"schema_version", "description", "notes", "flow_paths", "regulators", "side_paths", "not_modelled"}
# Inputs a tank's or bottle's own note is about (it would only repeat on the others).
NOTE_FIELDS = {"collapse_factor", "liquid_mass_kg", "outlet_cda_mm2", "minimum_temperature_k"}

ASSUMPTIONS = [
    ("Liquids", "LOX and RP-1 are incompressible, with properties (ThermoProp) at the tank set pressure and loaded "
                "temperature, except fuel downstream of regen. No boil-off during the burn."),
    ("Regen heating", "Fuel leaves the jacket at inlet + the configured rise (now 293 K to 375 K at every throttle level "
                      "and from the first instant, no thermal lag). The jacket's own dP uses its inlet properties, the "
                      "injector the heated ones. Combustion uses tank-temperature fuel: the regen heat comes from the "
                      "chamber gas, so the map is not rebuilt with hot fuel."),
    ("Ullage collapse", "Collapse factor per tank: the COPV supplies the full pressurant flow, the ullage keeps "
                        "1 / collapse_factor of it. LOX 2.0, fuel 1.1: assumed values; Luca to supply "
                        "conservative and Vespula-heritage numbers."),
    ("Ullage gas", "Real-gas nitrogen (ThermoProp), solved in pressure-enthalpy. No heat transfer with the tank "
                   "walls; only the collapse factor removes gas."),
    ("COPV", "Real-gas nitrogen blowdown with no wall heat transfer (adiabatic), starting at 293 K. The identical bottles "
             "are one lumped volume (pressurization.copv_model; \"header\" models each bottle feeding a common header)."),
    ("Tanks", "Hydrostatic head and vehicle acceleration are neglected (as in vehicle_sim.py)."),
    ("Pressurant lines", "Each tube is a lumped equivalent CdA, A / sqrt(1 + f L / D) with the fully rough friction "
                         "factor, in series with the leg's valves and orifices (FullFlow SeriesCdA). No gas volume "
                         "or acoustic dynamics; adiabatic."),
    ("Gas restrictions", "FullFlow CompressibleOrifice (choked / subcritical) with the upstream real-gas R and gamma."),
    ("Liquid lines", "FullFlow DarcyWeisbach with flow inertia and Churchill friction; junctions have no volume."),
    ("Throttle valves", "OTV / FTV are ideal flow controllers: their CdA is solved to pass the thrust-target flow, "
                        "capped at their 50 mm2 max CdA. They move between profile points at the GNC 180 lbf/s slew limit. "
                        "Physical valve travel is available (open_time_s / close_time_s) and is off."),
    ("Engine", "Equilibrium combustion map in (Pc, MR) times c* and Cf efficiencies that ramp linearly with Pc "
               "(c* 0.80 at 100 psia to 0.85 at 300 psia, Cf 0.94 at 100 to 0.97 at 250, flat outside; Luca's "
               "parameters, clamped instead of his max(1, ...)). In the transient each efficiency is taken at the "
               "previous step's Pc (one 0.05 s lag, under 0.1 % effect). Target MR 1.25. The chamber has no volume; propellant ignites as soon as it arrives."),
    ("Torch", "The current torch.network declares fuel FIMV -> FICK -> FIORF -> injector, oxidizer "
              "OIMV -> OICK -> OIORF -> OTHX -> injector, a shared torch chamber and throat, and purge "
              "TFPV -> FIPCK -> FIPORF into FIORF Outlet. OTHX is a tube with prescribed dT versus time "
              "since torch opening. The independent Hopper_Torch_Purge.py solves cold ideal-gas purge "
              "with a finite chamber and bidirectional geometric throat, no combustion or liquid displacement. "
              "Full vehicle coupling remains pending; keep torch.enabled false. Older "
              "configs with torch.branches retain the legacy direct-main-chamber approximation."),
    ("Start / shutdown", "GUESS sequence: OTV/FTV start opening at 0.30 s and take 0.20 s; the same ramp closes them "
                         "after the last profile point; 0.30 s of tail-off is simulated. The solver takes 0.01 s "
                         "steps through both ramps."),
    ("Bang-bang", "Each tank's solenoid closes above set + upper and opens below set - lower, sampled every "
                  "0.05 s, so pressure overshoots the band slightly. Solenoids switch instantly."),
    ("Regulator", "PR-REG is ideal: outlet = min(COPV, ambient + 750 psig) at the COPV gas temperature, no droop or lock-up."),
    ("Purges", "Purge GN2 is drawn from the COPV; OPV joins the LOX line below OTV, FPV the fuel line below FTV, "
               "Torch purge TFPV tees into FIORF Outlet and is solved only in the independent cold-gas case. Main purge "
               "gas only sees its feed-line back pressure; it is not added to the "
               "line or chamber flow. All purges are closed during the nominal burn."),
    ("Relief valves", "Opening rises linearly from crack to full-open dP (default 110 % of crack)."),
    ("Side paths", "Vent, relief, purge and drain flows only leave their source (no back-flow)."),
    ("Vent time", "Adiabatic blowdown of the ullage gas through the vent's CdA from set + upper limit to "
                  "1 psig, for a full tank (loaded ullage) and an empty tank. No LOX boil-off during the vent."),
    ("LOX drain", "OFD (same Jaksa valve as OVNT) plus any lox_drain.extra_elements as one CdA; ullage held at "
                  "each test pressure; outlet pressure = ullage + liquid head over a uniform cross-section "
                  "(tank volume / height). Incompressible LOX, no boil-off or flashing; the vented case is "
                  "optimistic because near-saturated LOX will partly flash in the valve."),
    ("Ambient", "Constant 14.8 psia (no altitude change during a hop)."),
    ("Hardware not modelled", "Quick disconnects (disconnected in flight), oxygen diffuser, torch spark "
                              "igniter, sensors."),
]

OPEN_ITEMS = [
    ("OBB / FBB CdA", "Model uses 18.68 mm2 for both. Master sheet: OBB = 2x PLV34 = 2.14 mm2, FBB = 1x PLV34 = "
                      "1.07 mm2. Which is right?"),
    ("Tank volumes", "Updated 2026-10-30 to the tank dimensions table: 15.867 L LOX / 17.929 L fuel, "
                     "16.00 / 13.00 kg loads. The master sheet's 19.88 / 15.18 L are out of date."),
    ("Line lengths and pressurant tube", "From Vedanth's 2026-11-05 config: LOX 48 + 10 in, fuel 12 + 12 in, pressurant "
                                         "1/4 in x 0.035 in, 48 in. Earlier: LOX 1.0 + 0.2032 m, fuel 0.3048 + 0.2032 m, "
                                         "pressurant 4.92 mm ID, 1.0 m. Which are measured?"),
    ("Regen outlet temperature", "Model uses 375 K (user, 2026-11-04). Vedanth's config had +51 K at 150 psia, +44 K at "
                                 "300 psia (344 / 337 K). Which is right? Should combustion see the hot fuel?"),
    ("ORV sizing", "RegO PRV9432T550 has a 0.062 in orifice (about 1.95 mm2 geometric); the sheet requires "
                   "30.8 mm2 CdA (20.5 mm2 minimum). The chosen part looks undersized; check the RegO flow table."),
    ("FRV set pressure", "Sheet FRV tab says 525 psi; the part (PRV9432T550) is 550 psig. Model uses 550."),
    ("OVNT vent time", "With the Jaksa 9.4 mm2, an empty 15.867 L LOX tank takes about 51 s to vent from 505 psia to "
                       "1 psig (requirement 60 s)."),
    ("Start-up tank spike", "Both tanks start at set pressure with the bang-bang solenoid open, so the first step "
                            "overshoots to about 570 psia, past the 565 psia relief crack. Start with "
                            "bang_bang_initial_state = 0?"),
    ("Vendor flow data", "No published Cv found for Marotta PLV34, Swagelok SS-4R3A, Check-All BUDSSTF.125SS or the "
                         "UNOX 3/8 NO solenoid. Need datasheets or flow tests."),
    ("LOX drain path", "Drain time uses only the OFD valve (9.4 mm2). OQD and the GSE drain line will slow it; add "
                       "them to lox_drain.extra_elements. Drain exit height below the tank and the drain ullage "
                       "pressure are GUESSes."),
    ("Nozzle note", "LinearizedNozzle docstring says Pc = Pa + 1 psi; nozzle_linear_pressure_drop_psi is 0.5."),
    ("Collapse factor", "LOX 2.0 / fuel 1.1 are assumed. Waiting on Luca for conservative and Vespula values."),
]


def _status(text):
    upper = text.upper()
    return "GUESS" if "GUESS" in upper else "EDIT" if "EDIT" in upper else ""


def _format(value):
    if isinstance(value, float):
        return float(f"{value:.6g}")
    if isinstance(value, (int, str, bool)):
        return value
    return "" if value is None else str(value)


def _scales():
    """{JSON field: factor to SI} from the same tables the loader uses."""
    scales = {key: scale for fields in GROUPS.values() for key, (_, scale) in fields.items()}
    scales.update({key: scale for key, scale in UNITS.values()})
    return scales


def _leaves(value, path, note=""):
    """(JSON path, value, note) for every input under a section."""
    if isinstance(value, dict):
        note = value.get("note", "") or note
        for key, item in value.items():
            if key not in ("note", "type"):
                yield from _leaves(item, f"{path}.{key}", note)
    elif isinstance(value, list) and value and all(isinstance(item, dict) and "name" in item for item in value):
        if path.startswith("bottles"):
            for item in value:
                yield from _leaves(item, f"{path}.{item['name']}")
        # Other lists of components (torch branches, drain extras) are on the Components sheet.
    else:
        yield path, value, note


def input_rows():
    """[section, name, as written, value (SI), note, status, owner] for every input in
    hopper_config.json. Flow paths and side paths are on the Components sheet."""
    notes, scales, rows = D.CONFIG.get("notes", {}), _scales(), []
    for section, content in D.CONFIG.items():
        if section in SECTIONS_ON_OTHER_SHEETS:
            continue
        for path, value, inherited in _leaves(content, section):
            key = path.rsplit(".", 1)[-1]
            if key == "time_s_thrust_n":
                continue
            note = notes.get(path, inherited if key in NOTE_FIELDS or section == "monte_carlo" else "")
            scale = scales.get(key, 1)
            if isinstance(value, bool) or not isinstance(value, (int, float)):
                si = _format(value)
            else:
                si = _format(float(value) * scale)
            rows.append([section, path, _format(value), si, note, _status(note), OWNERS.get(section, "")])
    return rows


def _sheet(workbook, title, header, rows, widths):
    sheet = workbook.create_sheet(title)
    sheet.append(header)
    for row in rows:
        sheet.append(list(row))
    bold, fill = Font(bold=True, color="FFFFFF"), PatternFill("solid", fgColor="305496")
    for cell in sheet[1]:
        cell.font, cell.fill = bold, fill
    guess = PatternFill("solid", fgColor="FFF2CC")
    status_column = header.index("Status") if "Status" in header else None
    for row in sheet.iter_rows(min_row=2):
        for cell in row:
            cell.alignment = Alignment(wrap_text=True, vertical="top")
        if status_column is not None and row[status_column].value == "GUESS":
            for cell in row:
                cell.fill = guess
    for k, width in enumerate(widths, start=1):
        sheet.column_dimensions[get_column_letter(k)].width = width
    sheet.freeze_panes = "A2"
    sheet.auto_filter.ref = sheet.dimensions
    return sheet


def write_workbook(directory=D.PROP_RESULTS_DIRECTORY):
    os.makedirs(directory, exist_ok=True)
    workbook = Workbook()
    workbook.remove(workbook.active)
    today = datetime.date.today().isoformat()

    inputs = [row + ["", ""] for row in input_rows()]
    _sheet(workbook, "Inputs",
           ["Section", "Name", "As written (units in the name)", "Value (SI)", "Note / source", "Status", "Owner team",
            "Owner", "Last confirmed"], inputs, [18, 44, 26, 14, 60, 9, 16, 14, 14])

    components = []
    for row in hardware_rows():
        status = _status(row[-1])
        path = row[0].lower()
        team = "Injectors" if row[2] == "Injector" else "Engines" if "regen" in row[1].lower() else "Fluids"
        components.append(row + [status, team if path != "not in the flow model" else "", "", ""])
    _sheet(workbook, "Components", COMPONENT_COLUMNS + ["Status", "Owner team", "Owner", "Last confirmed"],
           components, [40, 22, 18, 10, 46, 52, 9, 12, 14, 14])

    _sheet(workbook, "Assumptions", ["Topic", "Assumption"], ASSUMPTIONS, [22, 120])
    guesses = [(row[1], f"{row[4]} (value {row[3]})") for row in inputs if row[5] == "GUESS"]
    guesses += [(row[1], f"{row[5]} (CdA {row[3]} mm2)") for row in hardware_rows() if _status(row[-1]) == "GUESS"]
    _sheet(workbook, "Open items", ["Item", "Detail"],
           OPEN_ITEMS + [(f"GUESS: {name}", detail) for name, detail in guesses], [30, 120])

    about = workbook.create_sheet("About", 0)
    for line in [f"Hopper FullFlow parameters, generated {today} from hopper_config.json by Backend/Hopper_Parameters.py.",
                 "hopper_config.json is the source of truth: edit values there, then re-run this script.",
                 "Status GUESS (yellow rows) = not confirmed by the owning team. EDIT = expected to change.",
                 "Sheets: Inputs, Components (every P&ID item), Assumptions, Open items."]:
        about.append([line])
    about.column_dimensions["A"].width = 110

    path = os.path.join(directory, "Hopper_Parameters.xlsx")
    workbook.save(path)
    return path


if __name__ == "__main__":
    print("wrote", write_workbook())
