import fs from "node:fs/promises";
import path from "node:path";
import { fileURLToPath } from "node:url";

const artifactPath = "file:///C:/Program%20Files/WindowsApps/OpenAI.CodexPrimaryRuntime.v26-1112-737-0_26.1112.737.0_x64__3k8sg7r9htsxt/dependencies/node/node_modules/@oai/artifact-tool/dist/artifact_tool.mjs";
const { FileBlob, SpreadsheetFile, Workbook } = await import(artifactPath);

const here = path.dirname(fileURLToPath(import.meta.url));
const repoRoot = path.resolve(here, "../../../..");
const outputDir = path.join(here, "outputs", "valve_sizing_outputs");
const masterPath = path.join(repoRoot, "sims", "hopper_fullflow", "Outputs", "Hopper_Parameters.xlsx");
const config = JSON.parse(await fs.readFile(path.join(here, "valve_sizing_config.json"), "utf8"));

function parseCsv(text) {
  const lines = text.trim().split(/\r?\n/);
  const parseLine = (line) => {
    const result = []; let value = ""; let quote = false;
    for (let i = 0; i < line.length; i += 1) {
      const c = line[i];
      if (c === '"') { if (quote && line[i + 1] === '"') { value += '"'; i += 1; } else quote = !quote; }
      else if (c === "," && !quote) { result.push(value); value = ""; }
      else value += c;
    }
    result.push(value); return result;
  };
  const headers = parseLine(lines[0]);
  return lines.slice(1).filter(Boolean).map((line) => Object.fromEntries(headers.map((h, i) => [h, parseLine(line)[i] ?? ""])));
}
function number(v) { const n = Number(v); return Number.isFinite(n) ? n : v; }
async function csv(name) { return parseCsv(await fs.readFile(path.join(outputDir, name), "utf8")); }
const ventSummary = await csv("vent_summary.csv");
const ventSweep = await csv("vent_cda_sweep.csv");
const ventHistory = await csv("vent_history.csv");
const drainSummary = await csv("ofd_drain_summary.csv");
const drainSweep = await csv("ofd_drain_cda_sweep.csv");
const fillSummary = await csv("ofd_fill_summary.csv");
const fillSweep = await csv("ofd_fill_cda_sweep.csv");

const navy = "#1F4E78", blue = "#D9EAF7", pale = "#F3F7FA", yellow = "#FFF2CC", text = "#1F2937", border = "#B8C6D1";
const titleFmt = { font: { name: "Arial", size: 15, bold: true, color: text } };
const sectionFmt = { fill: navy, font: { name: "Arial", size: 10, bold: true, color: "#FFFFFF" }, horizontalAlignment: "left", verticalAlignment: "center" };
const headerFmt = { fill: "#305E8F", font: { name: "Arial", size: 10, bold: true, color: "#FFFFFF" }, horizontalAlignment: "center", verticalAlignment: "center", wrapText: true };
const bodyFmt = { font: { name: "Arial", size: 10, color: text }, verticalAlignment: "center" };
function styleTable(sheet, range, headerRows = 1) {
  sheet.getRange(range).format.borders = { preset: "all", style: "thin", color: border };
  const [start, end] = range.split(":");
  const row = start.match(/\d+/)[0];
  const firstCol = start.match(/[A-Z]+/)[0], lastCol = end.match(/[A-Z]+/)[0];
  sheet.getRange(`${firstCol}${row}:${lastCol}${row}`).format = headerFmt;
  sheet.getRange(range).format.verticalAlignment = "center";
}
function addTitle(sheet, title, subtitle = "") {
  sheet.showGridLines = false;
  sheet.getRange("A2:H2").merge(); sheet.getRange("A2").values = [[title]]; sheet.getRange("A2").format = titleFmt;
  if (subtitle) { sheet.getRange("A3:H3").merge(); sheet.getRange("A3").values = [[subtitle]]; sheet.getRange("A3").format = { font: { name: "Arial", size: 10, italic: true, color: "#55616D" } }; }
}
function addSection(sheet, cell, title, endCol = "H") {
  const row = Number(cell.match(/\d+/)[0]); const col = cell.match(/[A-Z]+/)[0];
  sheet.getRange(`${col}${row}:${endCol}${row}`).merge(); sheet.getRange(`${col}${row}`).values = [[title]]; sheet.getRange(`${col}${row}`).format = sectionFmt;
}
function writeKeyValue(sheet, startRow, rows) {
  sheet.getRange(`A${startRow}:D${startRow + rows.length - 1}`).values = rows.map(r => [r[0], r[1], r[2] ?? "", r[3] ?? ""]);
  sheet.getRange(`A${startRow}:D${startRow + rows.length - 1}`).format = bodyFmt;
  sheet.getRange(`A${startRow}:D${startRow + rows.length - 1}`).format.borders = { preset: "all", style: "thin", color: border };
  sheet.getRange(`A${startRow}:A${startRow + rows.length - 1}`).format.font = { name: "Arial", size: 10, bold: true, color: text };
}
function toRows(data, fields) { return data.map(d => fields.map(f => number(d[f]))); }

// New standalone sizing workbook.
const wb = Workbook.create();
const summary = wb.worksheets.add("Summary");
const inputs = wb.worksheets.add("Inputs");
const results = wb.worksheets.add("Results");
const pressure = wb.worksheets.add("Pressure history");
const assumptions = wb.worksheets.add("Assumptions");

addTitle(summary, "Hopper valve sizing", "Simulation output for OVNT, FVNT, and OFD. Values are trade-study outputs, not purchase approvals.");
addSection(summary, "A5", "Requirements and screening criteria", "D");
writeKeyValue(summary, 6, [
  ["Minimum component pressure rating", 625, "psia", "500 psia tank MEOP × 1.25 rating factor"],
  ["Minimum temperature rating", 80, "K or colder", "LOX saturation temperature with 10 K margin"],
  ["OVNT / FVNT vent endpoint", 1, "psig", "Target tank pressure after vent"],
  ["OFD drain endpoint", 1, "% liquid residual", "Drain screening endpoint"],
  ["OFD fill endpoint", 90, "% liquid fill", "Trade-study endpoint"],
]);
addSection(summary, "A13", "Results", "H");
summary.getRange("A14:H17").values = [
  ["Component", "Function", "CdA [mm²]", "Result", "Result units", "Basis", "Decision use", "Source"],
  ["OVNT", "LOX ullage vent", 9.4, number(ventSummary[0].vent_time_to_1psig_s), "s", "EOL blowdown", "Compare valve size vs vent time", "Vent results"],
  ["FVNT", "Fuel ullage vent", 15.0, number(ventSummary[1].vent_time_to_1psig_s), "s", "EOL blowdown", "Compare valve size vs vent time", "Vent results"],
  ["OFD", "LOX fill", number(fillSummary[0].example_OFD_CdA_mm2), number(fillSummary[0].liquid_transfer_time_s), "s liquid transfer", "250 psig LFS, fixed restrictions", "Compare valve size vs fill time", "OFD results"],
];
styleTable(summary, "A14:H17"); summary.getRange("C15:E17").format.numberFormat = "0.0";
summary.getRange("A20:H20").merge(); summary.getRange("A20").values = [["Use the pressure history tab for the pressure-versus-time graph. Inputs and assumptions are traceable in the adjacent tabs."]]; summary.getRange("A20").format = { font: { name: "Arial", size: 10, italic: true, color: "#55616D" } };
summary.getRange("A:A").format.columnWidth = 27; summary.getRange("B:B").format.columnWidth = 24; summary.getRange("C:E").format.columnWidth = 16; summary.getRange("F:H").format.columnWidth = 25;

addTitle(inputs, "Valve sizing inputs", "Editable settings belong in valve_sizing_config.json; this sheet records the values used for this simulation run.");
inputs.getRange("A5:D5").values = [["Input", "Value", "Units", "Source / basis"]]; styleTable(inputs, "A5:D5");
const inputRows = [
  ["Time step", config.time_step_s, "s", "valve_sizing_config.json"],
  ["Vent target pressure", config.vent_end_gauge_psig, "psig", "valve_sizing_config.json"],
  ["Vent target time", config.vent_target_time_s, "s", "valve_sizing_config.json"],
  ["OFD drain target time", config.drain_target_time_s, "s", "valve_sizing_config.json"],
  ["LFS supply pressure", config.lfs_supply_gauge_pressure_psig, "psig", "LOX Fill GFSSP Sim input"],
  ["LFS path CdA", config.ground_fill_path_cda_mm2, "mm²", "LOX Fill GFSSP Sim input"],
  ["OQD CdA", config.oqd_cda_mm2, "mm²", "Configurable provisional estimate"],
  ["Fill tank ullage pressure", config.fill_tank_ullage_gauge_pressure_psig, "psig", "valve_sizing_config.json"],
  ["Fill start fraction", config.fill_start_liquid_fraction, "fraction", "valve_sizing_config.json"],
  ["Fill end fraction", config.fill_end_liquid_fraction, "fraction", "valve_sizing_config.json"],
  ["Chilldown / boiloff delay", config.chilldown_and_boiloff_delay_s, "s", "Provisional first-pass setting"],
  ["Component pressure rating factor", config.minimum_component_pressure_rating_factor, "factor", "valve_sizing_config.json"],
  ["Minimum temperature margin", config.minimum_temperature_margin_k, "K", "valve_sizing_config.json"],
  ["LOX density", config.liquid_properties.LOX_density_kg_m3, "kg/m³", "valve_sizing_config.json"],
  ["RP-1 density", config.liquid_properties.RP1_density_kg_m3, "kg/m³", "valve_sizing_config.json"],
];
writeKeyValue(inputs, 6, inputRows); inputs.getRange("B6:B20").format.fill = yellow;
inputs.getRange("A:A").format.columnWidth = 32; inputs.getRange("B:B").format.columnWidth = 16; inputs.getRange("C:C").format.columnWidth = 16; inputs.getRange("D:D").format.columnWidth = 34;
inputs.freezePanes.freezeRows(5);

addTitle(results, "Valve sizing results", "Vent results use the FullFlow EOL tank state. OFD fill and drain results are separate liquid-flow trade studies.");
addSection(results, "A5", "Vent results", "N");
const ventFields = ["valve", "fluid", "tank", "configured_CdA_mm2", "water_equivalent_Cv", "EOL_liquid_mass_kg", "EOL_liquid_volume_L", "EOL_ullage_volume_L", "EOL_ullage_fraction", "start_pressure_psia", "ullage_temperature_K", "vent_time_to_1psig_s", "CdA_required_for_target_mm2", "meets_60s_target"];
const ventHeaders = ["Valve", "Fluid", "Tank", "CdA [mm²]", "Cv", "EOL liquid [kg]", "EOL liquid [L]", "EOL ullage [L]", "EOL ullage", "Start P [psia]", "Ullage T [K]", "Vent time [s]", "CdA for 60 s [mm²]", "≤60 s"];
results.getRange("A6:N8").values = [ventHeaders, ...toRows(ventSummary.slice(0,2), ventFields)]; styleTable(results, "A6:N8"); results.getRange("D7:M8").format.numberFormat = "0.000";
addSection(results, "A11", "OFD drain results", "I");
const drainFields = ["valve", "mode", "configured_CdA_mm2", "water_equivalent_Cv", "ullage_gauge_pressure_psig", "LOX_load_kg", "drain_time_to_1pct_residual_s", "CdA_required_for_target_mm2", "meets_60s_target"];
const drainHeaders = ["Valve", "Mode", "CdA [mm²]", "Cv", "Ullage P [psig]", "LOX load [kg]", "Drain time [s]", "CdA for 60 s [mm²]", "≤60 s"];
results.getRange(`A12:I${12 + drainSummary.length}`).values = [drainHeaders, ...toRows(drainSummary, drainFields)]; styleTable(results, `A12:I${12 + drainSummary.length}`); results.getRange(`C13:H${12 + drainSummary.length}`).format.numberFormat = "0.000";
addSection(results, "A20", "OFD fill results", "N");
const fillFields = ["valve", "mode", "example_OFD_CdA_mm2", "LFS_path_CdA_mm2", "OQD_CdA_mm2", "combined_path_CdA_mm2", "LFS_supply_gauge_pressure_psig", "tank_ullage_gauge_pressure_psig", "LOX_mass_added_kg", "fill_start_fraction", "fill_end_fraction", "chilldown_and_boiloff_delay_s", "liquid_transfer_time_s", "total_fill_operation_time_s"];
const fillHeaders = ["Valve", "Mode", "OFD CdA [mm²]", "LFS CdA [mm²]", "OQD CdA [mm²]", "Combined CdA [mm²]", "Supply P [psig]", "Tank P [psig]", "LOX added [kg]", "Start fraction", "End fraction", "Delay [s]", "Transfer [s]", "Total [s]"];
results.getRange("A21:N22").values = [fillHeaders, ...toRows(fillSummary, fillFields)]; styleTable(results, "A21:N22"); results.getRange("C22:N22").format.numberFormat = "0.000";
addSection(results, "A25", "CdA trade-study data", "F");
results.getRange(`A26:C${26 + ventSweep.length}`).values = [["Vent valve", "CdA [mm²]", "Vent time [s]"], ...ventSweep.map(r => [r.valve, number(r.cda_mm2), number(r.vent_time_s)])]; styleTable(results, `A26:C${26 + ventSweep.length}`);
results.getRange(`E26:F${26 + fillSweep.length}`).values = [["OFD CdA [mm²]", "Liquid transfer time [s]"], ...fillSweep.map(r => [number(r.ofd_cda_mm2), number(r.liquid_transfer_time_s)])]; styleTable(results, `E26:F${26 + fillSweep.length}`);
results.getRange("A:N").format.columnWidth = 15; results.getRange("B:B").format.columnWidth = 25; results.freezePanes.freezeRows(5);

addTitle(pressure, "EOL tank pressure versus time", "EOL tank pressure from the nitrogen blowdown simulation. OVNT and FVNT use their configured CdA values.");
const reduced = [];
for (const valve of ["OVNT", "FVNT"]) {
  const rows = ventHistory.filter(r => r.valve === valve);
  const step = Math.max(1, Math.floor(rows.length / 40));
  for (let i = 0; i < rows.length; i += step) reduced.push([valve, number(rows[i].time_s), number(rows[i].tank_pressure_psia), number(rows[i].nitrogen_mdot_g_s)]);
  if (rows.length && (reduced.length === 0 || reduced[reduced.length - 1][1] !== number(rows[rows.length - 1].time_s))) { const r = rows[rows.length - 1]; reduced.push([valve, number(r.time_s), number(r.tank_pressure_psia), number(r.nitrogen_mdot_g_s)]); }
}
const ovnt = reduced.filter(r => r[0] === "OVNT"), fvnt = reduced.filter(r => r[0] === "FVNT");
const maxLen = Math.max(ovnt.length, fvnt.length);
const chartRows = [["Time [s]", "OVNT pressure [psia]", "FVNT pressure [psia]"]];
for (let i=0; i<maxLen; i+=1) chartRows.push([ovnt[i]?.[1] ?? null, ovnt[i]?.[2] ?? null, fvnt[i]?.[2] ?? null]);
pressure.getRange(`A5:C${5 + chartRows.length - 1}`).values = chartRows; styleTable(pressure, `A5:C${5 + chartRows.length - 1}`); pressure.getRange(`A6:C${5 + chartRows.length - 1}`).format.numberFormat = "0.00";
const chart = pressure.charts.add("line", pressure.getRange(`A5:C${5 + chartRows.length - 1}`));
chart.title = "EOL Tank Pressure vs Time"; chart.hasLegend = true; chart.legend.position = "top"; chart.xAxis = { title: { text: "Time [s]" } }; chart.yAxis = { title: { text: "Tank pressure [psia]" }, numberFormatCode: "0", numberFormatSourceLinked: false }; chart.setPosition("E5", "N27");
pressure.getRange("E30:H31").merge(); pressure.getRange("E30").values = [["Raw pressure history: columns A-C are a reduced chart data series. The full-resolution simulation history remains in vent_history.csv."]]; pressure.getRange("E30").format = { font: { name: "Arial", size: 10, italic: true, color: "#55616D" }, wrapText: true };
pressure.getRange("A:C").format.columnWidth = 20; pressure.freezePanes.freezeRows(5);

addTitle(assumptions, "Model assumptions", "Assumptions retained from the sizing configuration used for this run.");
assumptions.getRange("A5:B5").values = [["Assumption", "Model treatment"]]; styleTable(assumptions, "A5:B5");
const assRows = Object.entries(config.assumptions).map(([k,v]) => [k.replace(/_/g, " "), v]);
assumptions.getRange(`A6:B${5 + assRows.length}`).values = assRows; assumptions.getRange(`A6:B${5 + assRows.length}`).format = bodyFmt; assumptions.getRange(`A5:B${5 + assRows.length}`).format.borders = { preset: "all", style: "thin", color: border }; assumptions.getRange(`B6:B${5 + assRows.length}`).format.wrapText = true;
assumptions.getRange("A:A").format.columnWidth = 30; assumptions.getRange("B:B").format.columnWidth = 105; assumptions.getRange(`A6:B${5 + assRows.length}`).format.autofitRows();

for (const sheet of [summary, inputs, results, pressure, assumptions]) sheet.tabColor = navy;
wb.recalculate();
await fs.mkdir(outputDir, { recursive: true });
const sizingOut = await SpreadsheetFile.exportXlsx(wb); await sizingOut.save(path.join(outputDir, "Valve_Sizing.xlsx"));

// Extend the existing component master without changing its source-driven component rows.
const master = await SpreadsheetFile.importXlsx(await FileBlob.load(masterPath));
let masterSheet;
try {
  masterSheet = master.worksheets.getItem("Valve sizing");
  masterSheet.deleteAllDrawings();
  masterSheet.getRange("A1:N400").unmerge();
  masterSheet.getRange("A1:N400").clear({ applyTo: "all" });
} catch { masterSheet = master.worksheets.add("Valve sizing"); }
masterSheet.showGridLines = false; masterSheet.tabColor = navy;
addTitle(masterSheet, "Valve sizing summary", "Linked local simulation output for OVNT, FVNT, and OFD. See the Valve_Sizing workbook for full tables and assumptions.");
addSection(masterSheet, "A5", "Requirements", "D");
writeKeyValue(masterSheet, 6, [
  ["Minimum pressure rating", 625, "psia", "500 psia MEOP × 1.25 factor"],
  ["Minimum temperature rating", 80, "K or colder", "LOX service with 10 K margin"],
  ["Material compatibility", "LOX / GN₂ compatible", "", "Select per vendor material documentation"],
  ["Connection / line limit", "1/2 in OD tubing", "", "Existing LOX line after QD"],
]);
addSection(masterSheet, "A12", "Inputs and assumptions", "D");
writeKeyValue(masterSheet, 13, [
  ["OVNT EOL ullage", number(ventSummary[0].EOL_ullage_fraction), "fraction", "Derived from EOL liquid residual, not 10% loaded ullage"],
  ["FVNT EOL ullage", number(ventSummary[1].EOL_ullage_fraction), "fraction", "Derived from EOL liquid residual, not 10% loaded ullage"],
  ["LFS supply pressure", config.lfs_supply_gauge_pressure_psig, "psig", "Fixed fill study input"],
  ["LFS path CdA", config.ground_fill_path_cda_mm2, "mm²", "Fixed fill study input"],
  ["OQD CdA", config.oqd_cda_mm2, "mm²", "Configurable provisional estimate"],
  ["Fill delay", config.chilldown_and_boiloff_delay_s, "s", "Provisional chilldown / boiloff delay"],
]);
masterSheet.getRange("B13:B14").format.numberFormat = "0.0%";
addSection(masterSheet, "A21", "Results", "H");
masterSheet.getRange("A22:H25").values = [
  ["Component", "Configured CdA [mm²]", "Metric", "Result", "Units", "Condition", "Use", "Workbook link"],
  ["OVNT", 9.4, "Vent to 1 psig", number(ventSummary[0].vent_time_to_1psig_s), "s", "EOL LOX tank state", "Vent time trade study", "=HYPERLINK(\"#'Valve sizing'!A1\",\"Open valve sizing\")"],
  ["FVNT", 15.0, "Vent to 1 psig", number(ventSummary[1].vent_time_to_1psig_s), "s", "EOL fuel tank state", "Vent time trade study", "=HYPERLINK(\"#'Valve sizing'!A1\",\"Open valve sizing\")"],
  ["OFD", number(fillSummary[0].example_OFD_CdA_mm2), "LOX liquid transfer", number(fillSummary[0].liquid_transfer_time_s), "s", "250 psig LFS / 90% fill", "Fill time trade study", "=HYPERLINK(\"#'Valve sizing'!A1\",\"Open valve sizing\")"],
];
styleTable(masterSheet, "A22:H25"); masterSheet.getRange("B23:D25").format.numberFormat = "0.000";
addSection(masterSheet, "A28", "Pressure versus time", "C");
masterSheet.getRange(`A29:C${29 + chartRows.length}`).values = chartRows; styleTable(masterSheet, `A29:C${29 + chartRows.length}`); masterSheet.getRange(`A30:C${29 + chartRows.length}`).format.numberFormat = "0.00";
const masterChart = masterSheet.charts.add("line", masterSheet.getRange(`A29:C${29 + chartRows.length}`));
masterChart.title = "EOL Tank Pressure vs Time"; masterChart.hasLegend = true; masterChart.legend.position = "top"; masterChart.xAxis = { title: { text: "Time [s]" } }; masterChart.yAxis = { title: { text: "Tank pressure [psia]" }, numberFormatCode: "0", numberFormatSourceLinked: false }; masterChart.setPosition("E28", "N50");
masterSheet.getRange("A:A").format.columnWidth = 32; masterSheet.getRange("B:D").format.columnWidth = 18; masterSheet.getRange("E:H").format.columnWidth = 24; masterSheet.freezePanes.freezeRows(5);

const components = master.worksheets.getItem("Components");
components.getRange("K1").values = [["Valve sizing"]]; components.getRange("K1").format = headerFmt; components.getRange("K2:K59").format = bodyFmt; components.getRange("K1:K59").format.borders = { preset: "all", style: "thin", color: border }; components.getRange("K:K").format.columnWidth = 21;
for (const row of [33, 35, 36]) components.getRange(`K${row}`).formulas = [["=HYPERLINK(\"#'Valve sizing'!A1\",\"Open valve sizing\")"]];

master.recalculate();
const masterOut = await SpreadsheetFile.exportXlsx(master); await masterOut.save(masterPath);

for (const [book, sheetName, range, outName] of [[wb, "Summary", "A1:H22", "Valve_Sizing_preview.png"], [wb, "Pressure history", "A1:N32", "Valve_Sizing_pressure_preview.png"], [master, "Valve sizing", "A1:N52", "Hopper_Parameters_valve_sizing_preview.png"], [master, "Components", "A30:K38", "Hopper_Parameters_components_preview.png"]]) {
  const image = await book.render({ sheetName, range, scale: 1, format: "png" });
  await fs.writeFile(path.join(outputDir, outName), new Uint8Array(await image.arrayBuffer()));
}

const check = await master.inspect({ kind: "table", range: "Valve sizing!A1:H25", include: "values,formulas", tableMaxRows: 30, tableMaxCols: 10 });
const errors = await master.inspect({ kind: "match", searchTerm: "#REF!|#DIV/0!|#VALUE!|#NAME\\?|#N/A|#NUM!|#NULL!|#SPILL!|#CALC!", options: { useRegex: true, maxResults: 100 }, summary: "final formula error scan" });
console.log(check.ndjson); console.log(errors.ndjson); console.log(JSON.stringify({ sizingWorkbook: path.join(outputDir, "Valve_Sizing.xlsx"), masterWorkbook: masterPath }, null, 2));
