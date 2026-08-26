## Responsive, draggable UI for the visualizer. It emits intent; it owns no flight state.
extends Control

signal telemetry_file_selected(path: String)
signal reference_file_selected(path: String)
signal tcp_requested
signal export_requested
signal settings_requested
signal playback_step_requested(amount: int)
signal playback_toggle_requested
signal speed_requested(multiplier: float)
signal scrub_requested(index: float)
signal camera_mode_requested(mode: int)
signal cloud_toggle_requested
signal vector_toggle_requested(kind: String, enabled: bool)
signal orbit_requested(relative_motion: Vector2)

const TelemetryGraphScript = preload("res://Scripts/telemetry_graph.gd")
const TelemetrySchema = preload("res://Scripts/telemetry_schema.gd")

var database: Node
var hud: Label
var status: Label
var timeline: HSlider
var graph_widgets: Array[Control] = []
var popup_graph_widgets: Array[Control] = []
var telemetry_panel: PanelContainer
var playback_panel: PanelContainer
var view_panel: PanelContainer
var force_panel: PanelContainer
var graphs_panel: PanelContainer
var graph_list: VBoxContainer
var file_dialog: FileDialog
var _file_action := "telemetry"

func configure(data_source: Node) -> void:
	database = data_source

func _ready() -> void:
	set_anchors_and_offsets_preset(Control.PRESET_FULL_RECT)
	_build()
	resized.connect(_layout_ui)
	call_deferred("reset_layout")

func update_readout(row: PackedFloat64Array, index: int, initial_fuel: float, initial_ox: float) -> void:
	timeline.set_value_no_signal(index)
	var velocity := Vector3(row[TelemetrySchema.Column.VN], row[TelemetrySchema.Column.VE], row[TelemetrySchema.Column.VD]).length()
	var mass := row[TelemetrySchema.Column.FUEL_MASS] + row[TelemetrySchema.Column.OX_MASS]
	var twr := row[TelemetrySchema.Column.THRUST] / maxf(mass * 9.80665, 1.0)
	var saturation := "SATURATING" if absf(row[TelemetrySchema.Column.CONTROL]) > 0.95 else "tracking"
	hud.text = "TELEMETRY  %s\nt = %.2f s   sample %d / %d\nAltitude %.1f m   speed %.1f m/s\nThrust %.0f N   mass %.1f kg   TWR %.2f\nTVC / RCS %+.2f  %s" % [database.source_label, row[TelemetrySchema.Column.TIME], index + 1, database.rows.size(), row[TelemetrySchema.Column.UP], velocity, row[TelemetrySchema.Column.THRUST], mass, twr, row[TelemetrySchema.Column.CONTROL], saturation]
	for graph in graph_widgets + popup_graph_widgets:
		if is_instance_valid(graph):
			graph.cursor = index
			graph.queue_redraw()

func set_sample_count(count: int) -> void:
	timeline.max_value = maxi(count - 1, 0)

func show_status(message: String) -> void:
	status.text = message

func reset_layout() -> void:
	for panel in [telemetry_panel, playback_panel, view_panel, force_panel, graphs_panel]:
		panel.set_meta("user_moved", false)
	_layout_ui()

func _build() -> void:
	telemetry_panel = _panel()
	var telemetry_stack := _stack(telemetry_panel, "TELEMETRY  •  drag this card")
	hud = _label("", 15)
	hud.add_theme_color_override("font_color", Color("e9f7ff"))
	telemetry_stack.add_child(hud)
	status = _label("Demo data loaded — choose CSV or start TCP :8080", 12)
	status.add_theme_color_override("font_color", Color("a9c7d9"))
	telemetry_stack.add_child(status)

	playback_panel = _panel()
	var controls := _stack(playback_panel, "PLAYBACK & DATA  •  drag this card")
	var sources := HFlowContainer.new()
	controls.add_child(sources)
	_button(sources, "Load CSV", _choose_telemetry)
	_button(sources, "TCP :8080", func(): tcp_requested.emit())
	_button(sources, "Export CSV", func(): export_requested.emit())
	_button(sources, "Reference CSV", _choose_reference)
	_button(sources, "Settings", func(): settings_requested.emit())
	_button(sources, "Reset layout", reset_layout)
	var playback := HFlowContainer.new()
	controls.add_child(playback)
	_button(playback, "◀ Step", func(): playback_step_requested.emit(-1))
	_button(playback, "Play / Pause", func(): playback_toggle_requested.emit())
	_button(playback, "Step ▶", func(): playback_step_requested.emit(1))
	var speed_select := OptionButton.new()
	for choice in ["0.25x", "0.5x", "1x", "1.5x", "2x"]:
		speed_select.add_item(choice)
	speed_select.select(2)
	speed_select.item_selected.connect(func(selected: int): speed_requested.emit([0.25, 0.5, 1.0, 1.5, 2.0][selected]))
	playback.add_child(speed_select)
	timeline = HSlider.new()
	timeline.custom_minimum_size.x = 130
	timeline.step = 1.0
	timeline.value_changed.connect(func(value: float): scrub_requested.emit(value))
	playback.add_child(timeline)

	view_panel = _panel()
	var camera_controls := _stack(view_panel, "CAMERA  •  drag this card")
	var camera_buttons := HFlowContainer.new()
	camera_controls.add_child(camera_buttons)
	_button(camera_buttons, "Fly-by", func(): camera_mode_requested.emit(0))
	_button(camera_buttons, "Orbit (right drag)", func(): camera_mode_requested.emit(1))
	_button(camera_buttons, "Engine landing", func(): camera_mode_requested.emit(2))
	_button(camera_buttons, "Cloud", func(): cloud_toggle_requested.emit())
	camera_controls.add_child(_label("Hold right mouse button and drag in the 3D view to orbit.", 12))

	graphs_panel = _panel()
	var analysis_stack := _stack(graphs_panel, "ANALYSIS  •  drag this card")
	var tabs := TabContainer.new()
	tabs.size_flags_vertical = Control.SIZE_EXPAND_FILL
	analysis_stack.add_child(tabs)
	var graphs := ScrollContainer.new()
	graphs.name = "Graphs"
	graphs.size_flags_vertical = Control.SIZE_EXPAND_FILL
	tabs.add_child(graphs)
	graph_list = VBoxContainer.new()
	graph_list.size_flags_horizontal = Control.SIZE_EXPAND_FILL
	graphs.add_child(graph_list)
	for column in range(1, TelemetrySchema.COLUMN_COUNT):
		var graph := TelemetryGraphScript.new()
		graph.database = database
		graph.column = column
		graph.title = TelemetrySchema.display_name(column)
		graph.pop_out_requested.connect(_open_graph_window)
		graph_list.add_child(graph)
		graph_widgets.append(graph)
	var pid := VBoxContainer.new()
	pid.name = "Fluid P&ID"
	pid.add_child(_label("Fluid system MVP", 18))
	pid.add_child(_label("Fuel valve: OPEN\nOx valve: OPEN\nPressure, temperature, and mass-flow telemetry can be added to the schema.", 14))
	tabs.add_child(pid)
	var body := VBoxContainer.new()
	body.name = "Body Frame"
	body.add_child(_label("CG & tank frame", 18))
	body.add_child(_label("Cyan: CG\nGold: CP\nBlue: fuel\nGreen: oxidizer\nOrange: thrust\nGreen: velocity\nPink: RCS", 14))
	tabs.add_child(body)

	force_panel = _panel()
	var vectors := _stack(force_panel, "VECTOR DISPLAY  •  drag this card")
	var vector_toggles := HFlowContainer.new()
	vectors.add_child(vector_toggles)
	for definition in [["Thrust", "thrust"], ["Velocity", "velocity"], ["RCS", "rcs"]]:
		var toggle := CheckButton.new()
		toggle.text = definition[0]
		toggle.button_pressed = true
		toggle.toggled.connect(func(enabled: bool): vector_toggle_requested.emit(definition[1], enabled))
		vector_toggles.add_child(toggle)

	file_dialog = FileDialog.new()
	file_dialog.file_mode = FileDialog.FILE_MODE_OPEN_FILE
	file_dialog.access = FileDialog.ACCESS_FILESYSTEM
	file_dialog.filters = PackedStringArray(["*.csv ; Telemetry CSV"])
	file_dialog.file_selected.connect(_file_selected)
	add_child(file_dialog)

func _choose_telemetry() -> void:
	_file_action = "telemetry"
	file_dialog.popup_centered_ratio(0.75)

func _choose_reference() -> void:
	_file_action = "reference"
	file_dialog.popup_centered_ratio(0.75)

func _file_selected(path: String) -> void:
	if _file_action == "telemetry":
		telemetry_file_selected.emit(path)
	else:
		reference_file_selected.emit(path)

func _open_graph_window(source: Control) -> void:
	var window := Window.new()
	window.title = "%s — Telemetry Graph" % source.title
	window.size = Vector2i(900, 560)
	window.min_size = Vector2i(520, 300)
	add_child(window)
	var graph := TelemetryGraphScript.new()
	graph.database = database
	graph.column = source.column
	graph.title = source.title
	graph.line_color = source.line_color
	graph.set_anchors_and_offsets_preset(Control.PRESET_FULL_RECT)
	window.add_child(graph)
	popup_graph_widgets.append(graph)
	window.close_requested.connect(func(): popup_graph_widgets.erase(graph); window.queue_free())
	window.popup_centered()

func _panel() -> PanelContainer:
	var panel := PanelContainer.new()
	panel.add_theme_stylebox_override("panel", _panel_style())
	add_child(panel)
	_enable_dragging(panel)
	return panel

func _stack(panel: PanelContainer, heading: String) -> VBoxContainer:
	var stack := VBoxContainer.new()
	stack.mouse_filter = Control.MOUSE_FILTER_IGNORE
	panel.add_child(stack)
	var handle := _label(heading, 11)
	handle.custom_minimum_size.y = 20.0
	handle.add_theme_color_override("font_color", Color("7fb5cf"))
	handle.mouse_filter = Control.MOUSE_FILTER_IGNORE
	stack.add_child(handle)
	return stack

func _button(parent: Container, label: String, action: Callable) -> void:
	var button := Button.new()
	button.text = label
	button.pressed.connect(action)
	parent.add_child(button)

func _label(text: String, font_size: int) -> Label:
	var label := Label.new()
	label.text = text
	label.add_theme_font_size_override("font_size", font_size)
	label.autowrap_mode = TextServer.AUTOWRAP_WORD_SMART
	return label

func _layout_ui() -> void:
	if size.x <= 1.0 or size.y <= 1.0:
		return
	var scale := clampf(minf(size.x / 1600.0, size.y / 900.0), 0.72, 1.20)
	var margin := 12.0 * scale
	var left_width := maxf(310.0 * scale, minf(510.0 * scale, size.x * 0.46))
	var graph_width := maxf(290.0 * scale, minf(430.0 * scale, size.x * 0.34))
	_apply_layout(telemetry_panel, Vector2(margin, margin), Vector2(left_width * 0.78, maxf(146.0, 156.0 * scale)))
	_apply_layout(playback_panel, Vector2(margin, telemetry_panel.position.y + telemetry_panel.size.y + margin), Vector2(left_width, maxf(150.0, 158.0 * scale)))
	_apply_layout(force_panel, Vector2(margin, playback_panel.position.y + playback_panel.size.y + margin), Vector2(minf(left_width, 350.0 * scale), 74.0))
	_apply_layout(view_panel, Vector2(margin, size.y - maxf(104.0, 110.0 * scale) - margin), Vector2(left_width * 0.82, maxf(104.0, 110.0 * scale)))
	_apply_layout(graphs_panel, Vector2(size.x - graph_width - margin, margin), Vector2(graph_width, maxf(330.0 * scale, size.y - margin * 2.0)))
	graph_list.custom_minimum_size.x = maxf(280.0, graphs_panel.size.x - 30.0)
	for graph in graph_widgets:
		graph.set_card_width(graphs_panel.size.x - 30.0)

func _apply_layout(panel: Control, default_position: Vector2, panel_size: Vector2) -> void:
	panel.size = panel_size
	if not panel.get_meta("user_moved", false):
		panel.position = default_position
	_clamp_panel(panel)

func _enable_dragging(panel: Control) -> void:
	var state := {"active": false, "offset": Vector2.ZERO}
	panel.gui_input.connect(func(event: InputEvent):
		if event is InputEventMouseButton:
			var button := event as InputEventMouseButton
			if button.button_index == MOUSE_BUTTON_LEFT:
				if button.pressed and button.position.y <= 26.0:
					state.active = true
					state.offset = button.global_position - panel.global_position
					panel.grab_click_focus()
				elif not button.pressed:
					state.active = false
		elif event is InputEventMouseMotion and state.active:
			var motion := event as InputEventMouseMotion
			panel.position = motion.global_position - state.offset
			panel.set_meta("user_moved", true)
			_clamp_panel(panel)
	)

func _clamp_panel(panel: Control) -> void:
	panel.position = Vector2(clampf(panel.position.x, 0.0, maxf(0.0, size.x - panel.size.x)), clampf(panel.position.y, 0.0, maxf(0.0, size.y - panel.size.y)))

func _input(event: InputEvent) -> void:
	if event is InputEventMouseMotion and Input.is_mouse_button_pressed(MOUSE_BUTTON_RIGHT):
		orbit_requested.emit((event as InputEventMouseMotion).relative)

func _panel_style() -> StyleBoxFlat:
	var style := StyleBoxFlat.new()
	style.bg_color = Color("0c1722e6")
	style.border_color = Color("3b617a")
	style.set_border_width_all(1)
	style.corner_radius_top_left = 8
	style.corner_radius_top_right = 8
	style.corner_radius_bottom_left = 8
	style.corner_radius_bottom_right = 8
	style.content_margin_left = 8
	style.content_margin_right = 8
	style.content_margin_top = 6
	style.content_margin_bottom = 6
	return style
