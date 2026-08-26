## Composition root: connects independent data, playback, world, and UI modules.
extends Node

const TelemetryDatabaseScript = preload("res://Scripts/telemetry_database.gd")
const PlaybackControllerScript = preload("res://Scripts/playback_controller.gd")
const FlightWorldScript = preload("res://Scripts/flight_world.gd")
const DashboardScript = preload("res://Scripts/visualizer_dashboard.gd")
const DemoTelemetry = preload("res://Scripts/demo_telemetry.gd")
const TelemetrySchema = preload("res://Scripts/telemetry_schema.gd")

var database: Node
var playback: Node
var flight_world: Node
var dashboard: Control
var initial_fuel := 1.0
var initial_ox := 1.0
var cloud_visible := true

func _ready() -> void:
	database = TelemetryDatabaseScript.new()
	playback = PlaybackControllerScript.new()
	flight_world = FlightWorldScript.new()
	dashboard = DashboardScript.new()
	add_child(database)
	add_child(playback)
	add_child(flight_world)
	var layer := CanvasLayer.new()
	add_child(layer)
	dashboard.configure(database)
	layer.add_child(dashboard)
	_connect_modules()
	database.rows = DemoTelemetry.create_rows()
	database.data_changed.emit()

func _connect_modules() -> void:
	database.data_changed.connect(_refresh_data)
	playback.sample_changed.connect(_display_sample)
	dashboard.telemetry_file_selected.connect(_load_telemetry)
	dashboard.reference_file_selected.connect(_load_reference)
	dashboard.tcp_requested.connect(_start_tcp)
	dashboard.export_requested.connect(_export_csv)
	dashboard.settings_requested.connect(_show_settings)
	dashboard.playback_step_requested.connect(playback.step)
	dashboard.playback_toggle_requested.connect(playback.toggle_playing)
	dashboard.speed_requested.connect(playback.set_speed)
	dashboard.scrub_requested.connect(playback.scrub_to)
	dashboard.camera_mode_requested.connect(flight_world.set_camera_mode)
	dashboard.cloud_toggle_requested.connect(_toggle_cloud)
	dashboard.vector_toggle_requested.connect(flight_world.set_vector_enabled)
	dashboard.orbit_requested.connect(flight_world.orbit)

func _refresh_data() -> void:
	if database.rows.is_empty():
		return
	initial_fuel = maxf(database.rows[0][TelemetrySchema.Column.FUEL_MASS], 0.001)
	initial_ox = maxf(database.rows[0][TelemetrySchema.Column.OX_MASS], 0.001)
	playback.set_sample_count(database.rows.size())
	flight_world.set_trace(database.rows)
	dashboard.set_sample_count(database.rows.size())
	playback.emit_current_sample()

func _display_sample(index: int) -> void:
	var row: PackedFloat64Array = database.sample_at(index)
	if not row.is_empty():
		flight_world.apply_sample(index, row, initial_fuel, initial_ox)
		dashboard.update_readout(row, index, initial_fuel, initial_ox)

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

func _show_settings() -> void:
	dashboard.show_status("Settings: landing ellipse 14 × 8 m; Monte Carlo cloud 60 points; TCP archive cap 1000 rows.")

func _toggle_cloud() -> void:
	cloud_visible = not cloud_visible
	flight_world.set_cloud_visible(cloud_visible)

func _process(delta: float) -> void:
	flight_world.update_camera(delta)
