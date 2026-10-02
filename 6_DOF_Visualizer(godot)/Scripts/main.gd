## Composition root: connects independent data, playback, world, and UI modules.
extends Node

const TelemetryDatabaseScript = preload("res://Scripts/telemetry_database.gd")
const PlaybackControllerScript = preload("res://Scripts/playback_controller.gd")
const FlightWorldScript = preload("res://Scripts/flight_world.gd")
const DashboardScript = preload("res://Scripts/visualizer_dashboard.gd")
const DemoTelemetry = preload("res://Scripts/demo_telemetry.gd")
const TelemetrySchema = preload("res://Scripts/telemetry_schema.gd")
const TelemetryInterpolator = preload("res://Scripts/telemetry_interpolator.gd")
const AdministrationPanelScene = preload("res://Scenes/administration_panel.tscn")
const PdfReportExporter = preload("res://Scripts/pdf_report_exporter.gd")

var database: Node
var playback: Node
var flight_world: Node
var dashboard: Control
var interface_layer: CanvasLayer
var initial_fuel := 1.0
var initial_ox := 1.0
var administration_panel: Control

func _ready() -> void:
	database = TelemetryDatabaseScript.new()
	playback = PlaybackControllerScript.new()
	flight_world = FlightWorldScript.new()
	dashboard = DashboardScript.new()
	add_child(database)
	add_child(playback)
	add_child(flight_world)
	interface_layer = CanvasLayer.new()
	add_child(interface_layer)
	dashboard.configure(database)
	interface_layer.add_child(dashboard)
	_connect_modules()
	database.rows = DemoTelemetry.create_rows()
	database.data_changed.emit()

func _connect_modules() -> void:
	database.data_changed.connect(_refresh_data)
	playback.frame_changed.connect(_display_frame)
	dashboard.telemetry_file_selected.connect(_load_telemetry)
	dashboard.reference_file_selected.connect(_load_reference)
	dashboard.tcp_requested.connect(_start_tcp)
	dashboard.export_requested.connect(_export_csv)
	dashboard.administration_requested.connect(_open_administration)
	dashboard.playback_step_requested.connect(playback.step)
	dashboard.playback_toggle_requested.connect(playback.toggle_playing)
	dashboard.speed_requested.connect(playback.set_speed)
	dashboard.loop_enabled_requested.connect(playback.set_loop_enabled)
	dashboard.scrub_requested.connect(playback.scrub_to)
	dashboard.camera_mode_requested.connect(flight_world.set_camera_mode)
	dashboard.landing_ellipse_toggle_requested.connect(flight_world.set_landing_ellipse_visible)
	dashboard.monte_carlo_toggle_requested.connect(flight_world.set_monte_carlo_visible)
	dashboard.vector_toggle_requested.connect(flight_world.set_vector_enabled)
	dashboard.orbit_requested.connect(flight_world.orbit)

func _refresh_data() -> void:
	if database.rows.is_empty():
		return
	initial_fuel = maxf(database.rows[0][TelemetrySchema.Column.FUEL_MASS], 0.001)
	initial_ox = maxf(database.rows[0][TelemetrySchema.Column.OX_MASS], 0.001)
	playback.set_rows(database.rows)
	flight_world.set_trace(database.rows)
	dashboard.set_sample_count(database.rows.size())
	playback.emit_current_frame()

func _display_frame(start_index: int, end_index: int, interpolation: float) -> void:
	var row: PackedFloat64Array = database.sample_at(start_index)
	if end_index != start_index:
		row = TelemetryInterpolator.interpolate(row, database.sample_at(end_index), interpolation)
	
	if not row.is_empty():
		flight_world.apply_sample(start_index, row, initial_fuel, initial_ox)
		dashboard.update_readout(row, start_index, initial_fuel, initial_ox)
		
		# --- UPDATE THE LANDING ELLIPSE DYNAMICALLY PER FRAME ---
		var landing_pred := _calculate_predicted_landing(row)
		flight_world.set_landing_ellipse_params(
			landing_pred.center_east,
			landing_pred.center_north,
			landing_pred.semi_major,
			landing_pred.semi_minor,
			landing_pred.rotation
		)
		# --------------------------------------------------------

func _load_telemetry(path: String) -> void:
	if database.load_csv(path):
		playback.reset()
		dashboard.show_status("Loaded %s (%d readings)" % [path.get_file(), database.rows.size()])

func _load_reference(path: String) -> void:
	var reference_rows: Array[PackedFloat64Array] = database.read_csv_rows(path)
	flight_world.set_reference_trace(reference_rows)
	dashboard.show_status("Reference trajectory loaded (%d points)" % reference_rows.size())

func _start_tcp() -> void:
	if database.start_tcp(8080):
		dashboard.show_status("Listening for newline-delimited 17-column CSV on TCP :8080")

func _export_csv() -> void:
	var path := "user://visualizer_export_%s.csv" % Time.get_datetime_string_from_system().replace(":", "-")
	if database.export_csv(path):
		dashboard.show_status("Exported telemetry to %s" % path)

func _open_administration() -> void:
	if is_instance_valid(administration_panel):
		return
	administration_panel = AdministrationPanelScene.instantiate()
	administration_panel.configure(database)
	administration_panel.back_requested.connect(_close_administration)
	administration_panel.pdf_export_requested.connect(_export_pdf_report)
	administration_panel.palette_requested.connect(flight_world.set_overlay_palette)
	administration_panel.landing_ellipse_toggled.connect(flight_world.set_landing_ellipse_visible)
	administration_panel.monte_carlo_toggled.connect(flight_world.set_monte_carlo_visible)
	dashboard.hide()
	interface_layer.add_child(administration_panel)

func _close_administration() -> void:
	if is_instance_valid(administration_panel):
		administration_panel.queue_free()
	administration_panel = null
	dashboard.show()

func _export_pdf_report(path: String) -> void:
	if PdfReportExporter.export_report(path, database):
		dashboard.show_status("Exported graph report to %s" % path)
		if is_instance_valid(administration_panel):
			administration_panel.show_status("Exported graph report to %s" % path)
	else:
		dashboard.show_status("Could not export graph report to %s" % path)
		if is_instance_valid(administration_panel):
			administration_panel.show_status("Could not export graph report to %s" % path)

func _process(delta: float) -> void:
	flight_world.update_camera(delta)
	
	## Computes the instantaneous predicted landing point and shrinks uncertainty as altitude decreases.
func _calculate_predicted_landing(row: PackedFloat64Array) -> Dictionary:
	var north := row[TelemetrySchema.Column.NORTH]
	var east := row[TelemetrySchema.Column.EAST]
	var up := maxf(row[TelemetrySchema.Column.UP], 0.0)
	var vn := row[TelemetrySchema.Column.VN]
	var ve := row[TelemetrySchema.Column.VE]
	var vd := row[TelemetrySchema.Column.VD] # Positive is downward
	
	var pred_north := north
	var pred_east := east
	
	# If descending and above ground, project forward linearly (or ballistically)
	if vd > 0.1 and up > 0.0:
		var time_to_go := up / vd
		pred_north += vn * time_to_go
		pred_east += ve * time_to_go
	
	# Uncertainty ellipse shrinks as altitude decreases (converging to accuracy at touchdown)
	var altitude_fraction := clampf(up / 1000.0, 0.05, 1.0)
	var semi_major := maxf(14.0 * altitude_fraction, 3.0)
	var semi_minor := maxf(8.0 * altitude_fraction, 2.0)
	
	return {
		"center_east": pred_east,
		"center_north": pred_north,
		"semi_major": semi_major,
		"semi_minor": semi_minor,
		"rotation": 0.0
	}
