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
const MainMenuScript = preload("res://Scripts/main_menu.gd")
const PresentationSettingsScript = preload("res://Scripts/presentation_settings.gd")
const VehicleProfileStoreScript = preload("res://Scripts/vehicle_profile_store.gd")
const SettingsPageScene = preload("res://Scenes/settings_page.tscn")
const VehicleSelectorScene = preload("res://Scenes/vehicle_selector.tscn")

var database: Node
var playback: Node
var flight_world: Node
var dashboard: Control
var interface_layer: CanvasLayer
var main_menu: Control
var initial_fuel := 1.0
var initial_ox := 1.0
var administration_panel: Control
var presentation_settings: Node
var vehicle_profiles: Node
var settings_page: Control
var vehicle_selector: Control

func _ready() -> void:
	database = TelemetryDatabaseScript.new()
	playback = PlaybackControllerScript.new()
	flight_world = FlightWorldScript.new()
	dashboard = DashboardScript.new()
	presentation_settings = PresentationSettingsScript.new()
	vehicle_profiles = VehicleProfileStoreScript.new()
	add_child(database)
	add_child(playback)
	add_child(flight_world)
	add_child(presentation_settings)
	add_child(vehicle_profiles)
	presentation_settings.load_from_disk()
	vehicle_profiles.load_from_disk()
	interface_layer = CanvasLayer.new()
	add_child(interface_layer)
	dashboard.configure(database, flight_world)
	interface_layer.add_child(dashboard)
	main_menu = MainMenuScript.new()
	main_menu.configure(vehicle_profiles, presentation_settings)
	main_menu.load_mission_requested.connect(_open_visualizer)
	main_menu.settings_requested.connect(_open_settings)
	main_menu.vehicle_selector_requested.connect(_open_vehicle_selector)
	interface_layer.add_child(main_menu)
	dashboard.hide()
	_connect_modules()
	database.rows = DemoTelemetry.create_rows()
	database.data_changed.emit()
	_apply_vehicle_profile(vehicle_profiles.get_selected_profile())
	_apply_presentation(presentation_settings.snapshot())

func _connect_modules() -> void:
	database.data_changed.connect(_refresh_data)
	playback.frame_changed.connect(_display_frame)
	dashboard.telemetry_file_selected.connect(_load_telemetry)
	dashboard.reference_file_selected.connect(_load_reference)
	dashboard.tcp_requested.connect(_start_tcp)
	dashboard.export_requested.connect(_export_csv)
	dashboard.administration_requested.connect(_open_administration)
	dashboard.main_menu_requested.connect(_return_to_main_menu)
	dashboard.playback_step_requested.connect(playback.step)
	dashboard.playback_toggle_requested.connect(playback.toggle_playing)
	dashboard.speed_requested.connect(playback.set_speed)
	dashboard.loop_enabled_requested.connect(playback.set_loop_enabled)
	dashboard.scrub_requested.connect(playback.scrub_to)
	dashboard.landing_ellipse_toggle_requested.connect(flight_world.set_landing_ellipse_visible)
	dashboard.monte_carlo_toggle_requested.connect(flight_world.set_monte_carlo_visible)
	dashboard.terrain_load_requested.connect(flight_world.load_terrain)
	dashboard.terrain_reset_requested.connect(flight_world.restore_default_terrain)
	dashboard.vector_toggle_requested.connect(flight_world.set_vector_enabled)
	dashboard.orbit_requested.connect(flight_world.orbit)
	flight_world.environment_status_changed.connect(dashboard.show_status)
	presentation_settings.changed.connect(_apply_presentation)
	vehicle_profiles.selected_profile_changed.connect(_apply_vehicle_profile)

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

func _open_visualizer() -> void:
	if is_instance_valid(administration_panel):
		_close_administration()
	main_menu.hide()
	dashboard.show()
	# The dashboard is created while the Main Menu is visible.  Always begin a
	# mission from its first telemetry sample when entering the workspace.
	playback.restart_playing()
	dashboard.show_status("Demo mission loaded — choose CSV or load a reference trajectory.")

func _open_settings() -> void:
	if is_instance_valid(settings_page):
		return
	settings_page = SettingsPageScene.instantiate()
	settings_page.configure(presentation_settings)
	settings_page.back_requested.connect(_close_settings)
	main_menu.hide()
	interface_layer.add_child(settings_page)

func _close_settings() -> void:
	if is_instance_valid(settings_page):
		settings_page.queue_free()
	settings_page = null
	main_menu.show()

func _open_vehicle_selector() -> void:
	if is_instance_valid(vehicle_selector):
		return
	vehicle_selector = VehicleSelectorScene.instantiate()
	vehicle_selector.configure(vehicle_profiles, presentation_settings)
	vehicle_selector.back_requested.connect(_close_vehicle_selector)
	main_menu.hide()
	interface_layer.add_child(vehicle_selector)

func _close_vehicle_selector() -> void:
	if is_instance_valid(vehicle_selector):
		vehicle_selector.queue_free()
	vehicle_selector = null
	main_menu.show()

func _return_to_main_menu() -> void:
	if is_instance_valid(administration_panel):
		administration_panel.queue_free()
		administration_panel = null
	# Do not let an invisible visualization continue consuming its playback
	# clock behind the Main Menu.  Rewinding also makes the next Load Mission
	# entry deterministic.
	playback.rewind_and_pause()
	dashboard.hide()
	main_menu.show()

func _apply_presentation(colors: Dictionary) -> void:
	if is_instance_valid(flight_world):
		flight_world.set_visual_colors(colors)
	if is_instance_valid(main_menu):
		main_menu.apply_presentation(colors)
	if is_instance_valid(dashboard):
		dashboard.apply_presentation(colors)

func _apply_vehicle_profile(profile: Dictionary) -> void:
	if not is_instance_valid(flight_world):
		return
	var model: Node3D = vehicle_profiles.instantiate_model(profile)
	if model != null:
		flight_world.set_vehicle_model(model)

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
	flight_world.update_view_cameras(delta, dashboard.get_view_cameras())
