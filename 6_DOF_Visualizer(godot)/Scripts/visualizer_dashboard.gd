## Fixed visualization workspace. This control owns presentation only; flight state remains in FlightWorld.
extends Control

signal telemetry_file_selected(path: String)
signal reference_file_selected(path: String)
signal tcp_requested
signal export_requested
signal administration_requested
signal main_menu_requested
signal playback_step_requested(amount: int)
signal playback_toggle_requested
signal speed_requested(multiplier: float)
signal loop_enabled_requested(enabled: bool)
signal scrub_requested(index: float)
signal landing_ellipse_toggle_requested(visible: bool)
signal monte_carlo_toggle_requested(visible: bool)
signal terrain_load_requested(latitude: float, longitude: float, api_key: String, area_size_km: float)
signal terrain_reset_requested
signal vector_toggle_requested(kind: String, enabled: bool)
signal orbit_requested(relative_motion: Vector2)

const TelemetryGraphScript = preload("res://Scripts/telemetry_graph.gd")
const TelemetrySchema = preload("res://Scripts/telemetry_schema.gd")
const VIEW_NAMES := {0: "FLY-BY", 1: "MOUSE CAMERA", 2: "ENGINE LANDING"}

var database: Node
var flight_world: Node
var hud: Label
var status: Label
var timeline: HSlider
var graph_widgets: Array[Control] = []
var popup_graph_widgets: Array[Control] = []
var top_bar: PanelContainer
var top_dropdown_band: PanelContainer
var bottom_band: PanelContainer
var telemetry_panel: PanelContainer
var playback_panel: PanelContainer
var view_panel: PanelContainer
var force_panel: PanelContainer
var reference_panel: PanelContainer
var monte_carlo_panel: PanelContainer
var graphs_panel: PanelContainer
var graphs_tab: Button
var graph_list: VBoxContainer
var collapsed_graph_list: VBoxContainer
var workspace: Control
var view_grid: Control
var file_dialog: FileDialog
var terrain_latitude_field: LineEdit
var terrain_longitude_field: LineEdit
var terrain_key_field: LineEdit
var terrain_area_selector: OptionButton
var _file_action := "telemetry"
var _graphs_open := false
var _active_view_modes: Array[int] = [0, 1, 2]
var _primary_view_mode := 0
var _view_buttons: Dictionary = {}
var _view_cameras: Dictionary = {}
var _view_panes: Dictionary = {}
var _viewports_ready := false
var _presentation_colors: Dictionary = {"ui_base": Color("062e49"), "ui_highlight": Color("66cbed")}

func configure(data_source: Node, world: Node) -> void:
	database = data_source
	flight_world = world

func _ready() -> void:
	set_anchors_and_offsets_preset(Control.PRESET_FULL_RECT)
	_build()
	resized.connect(_layout_ui)
	call_deferred("_create_viewports")
	call_deferred("_layout_ui")

func update_readout(row: PackedFloat64Array, index: int, initial_fuel: float, initial_ox: float) -> void:
	timeline.set_value_no_signal(index)
	var velocity := Vector3(row[TelemetrySchema.Column.VN], row[TelemetrySchema.Column.VE], row[TelemetrySchema.Column.VD]).length()
	var mass := row[TelemetrySchema.Column.FUEL_MASS] + row[TelemetrySchema.Column.OX_MASS]
	var twr := row[TelemetrySchema.Column.THRUST] / maxf(mass * 9.80665, 1.0)
	hud.text = "t %.2f s  |  %d / %d\nALT %.1f m   V %.1f m/s\nTHR %.0f N   TWR %.2f" % [row[TelemetrySchema.Column.TIME], index + 1, database.rows.size(), row[TelemetrySchema.Column.UP], velocity, row[TelemetrySchema.Column.THRUST], twr]
	for graph in graph_widgets + popup_graph_widgets:
		if is_instance_valid(graph):
			graph.cursor = index
			graph.queue_redraw()

func set_sample_count(count: int) -> void:
	timeline.max_value = maxi(count - 1, 0)

func show_status(message: String) -> void:
	status.text = message

func get_view_cameras() -> Dictionary:
	return _view_cameras

## Applies only palette values.  Geometry, margins, and the original dashboard
## component styles remain owned by this scene so Appearance cannot reflow it.
func apply_presentation(colors: Dictionary) -> void:
	_presentation_colors = colors.duplicate(true)
	for panel in [top_bar, top_dropdown_band, bottom_band, telemetry_panel, playback_panel, view_panel, force_panel, reference_panel, monte_carlo_panel, graphs_panel]:
		if is_instance_valid(panel):
			panel.add_theme_stylebox_override("panel", _panel_style(panel == telemetry_panel))
	if is_instance_valid(top_bar) and top_bar.get_child_count() > 0:
		for control in top_bar.get_child(0).get_children():
			if control is Button:
				_apply_navigation_style(control as Button)
	if is_instance_valid(graphs_tab):
		_apply_navigation_style(graphs_tab)
	for mode in _view_panes:
		var pane := _view_panes[mode] as PanelContainer
		if pane != null:
			pane.add_theme_stylebox_override("panel", _pane_style())
	_update_camera_button_styles()

func _build() -> void:
	workspace = Control.new()
	workspace.name = "CameraWorkspace"
	workspace.mouse_filter = Control.MOUSE_FILTER_IGNORE
	add_child(workspace)
	view_grid = Control.new()
	view_grid.mouse_filter = Control.MOUSE_FILTER_PASS
	workspace.add_child(view_grid)

	top_bar = _panel()
	var navigation := HBoxContainer.new()
	navigation.add_theme_constant_override("separation", 8)
	top_bar.add_child(navigation)
	_navigation_button(navigation, "MENU", func(): main_menu_requested.emit())
	_navigation_button(navigation, "FORCES", func(): _toggle_dropdown(force_panel))
	_navigation_button(navigation, "REFERENCE", func(): _toggle_dropdown(reference_panel))
	_navigation_button(navigation, "ANIMATIONS", Callable(), "Animation controls will be added after their telemetry inputs are defined.")
	var spacer := Control.new()
	spacer.size_flags_horizontal = Control.SIZE_EXPAND_FILL
	navigation.add_child(spacer)
	_navigation_button(navigation, "EXPORT", func(): administration_requested.emit())
	top_dropdown_band = _panel()
	top_dropdown_band.visible = false

	force_panel = _panel()
	force_panel.visible = false
	var forces := _stack(force_panel, "FORCES")
	var vector_toggles := HFlowContainer.new()
	_flow_layout(vector_toggles)
	forces.add_child(vector_toggles)
	for definition in [["Thrust", "thrust"], ["Velocity", "velocity"]]:
		var toggle := CheckButton.new()
		toggle.text = definition[0]
		toggle.button_pressed = true
		toggle.toggled.connect(func(enabled: bool): vector_toggle_requested.emit(definition[1], enabled))
		vector_toggles.add_child(toggle)
	_disabled_button(vector_toggles, "Drag", "Drag-force telemetry is not in the current data contract.")
	_disabled_button(vector_toggles, "Slosh", "Slosh-force telemetry is not in the current data contract.")

	reference_panel = _panel()
	reference_panel.visible = false
	var reference := _stack(reference_panel, "REFERENCE GEOMETRY")
	var reference_actions := HFlowContainer.new()
	_flow_layout(reference_actions)
	reference.add_child(reference_actions)
	_button(reference_actions, "Load reference CSV", _choose_reference)
	_disabled_button(reference_actions, "Landing ellipse", "Landing-ellipse settings remain in the Settings page.")
	_disabled_button(reference_actions, "Projection options", "Reference-layer controls will be added with the final settings plan.")
	var terrain_controls := HFlowContainer.new()
	_flow_layout(terrain_controls)
	reference.add_child(terrain_controls)
	terrain_latitude_field = _terrain_field("Latitude", false)
	terrain_longitude_field = _terrain_field("Longitude", false)
	terrain_key_field = _terrain_field("OpenTopo key (session only)", true)
	terrain_area_selector = OptionButton.new()
	terrain_area_selector.custom_minimum_size = Vector2(112.0, 30.0)
	terrain_area_selector.tooltip_text = "Width and depth of the square terrain area centered on the selected coordinates."
	terrain_area_selector.add_item("1 × 1 km", 1)
	terrain_area_selector.add_item("2 × 2 km", 2)
	terrain_area_selector.add_item("5 × 5 km", 5)
	terrain_controls.add_child(terrain_latitude_field)
	terrain_controls.add_child(terrain_longitude_field)
	terrain_controls.add_child(terrain_key_field)
	terrain_controls.add_child(terrain_area_selector)
	_button(terrain_controls, "Load terrain", _request_terrain)
	_button(terrain_controls, "Restore test range", func(): terrain_reset_requested.emit())
	reference.add_child(_label("Terrain uses only this session's coordinates and key; nothing is saved.", 10))

	telemetry_panel = _panel(true)
	var telemetry_stack := _stack(telemetry_panel, "FLIGHT STATUS")
	hud = _label("Awaiting telemetry", 12)
	hud.add_theme_color_override("font_color", Color("e9f7ff"))
	telemetry_stack.add_child(hud)
	status = _label("Demo mission loaded", 10)
	status.add_theme_color_override("font_color", Color("a9c7d9"))
	status.autowrap_mode = TextServer.AUTOWRAP_WORD_SMART
	telemetry_stack.add_child(status)

	bottom_band = _panel()
	playback_panel = _panel()
	var controls := _stack(playback_panel, "PLAYBACK & DATA")
	var sources := HFlowContainer.new()
	_flow_layout(sources)
	controls.add_child(sources)
	_button(sources, "Load CSV", _choose_telemetry)
	_button(sources, "TCP :8080", func(): tcp_requested.emit())
	_button(sources, "Reference", _choose_reference)
	_button(sources, "Export CSV", func(): export_requested.emit())
	var playback_controls := HFlowContainer.new()
	_flow_layout(playback_controls)
	controls.add_child(playback_controls)
	_button(playback_controls, "◀ Step", func(): playback_step_requested.emit(-1))
	_button(playback_controls, "Play / Pause", func(): playback_toggle_requested.emit())
	_button(playback_controls, "Step ▶", func(): playback_step_requested.emit(1))
	_disabled_button(playback_controls, "To Live", "Live-session navigation will be enabled with the finalized TCP session flow.")
	var speed_select := OptionButton.new()
	for choice in ["0.25x", "0.5x", "1x real-time", "1.5x", "2x", "4x", "8x"]:
		speed_select.add_item(choice)
	speed_select.select(2)
	speed_select.item_selected.connect(func(selected: int): speed_requested.emit([0.25, 0.5, 1.0, 1.5, 2.0, 4.0, 8.0][selected]))
	playback_controls.add_child(speed_select)
	var loop_toggle := CheckButton.new()
	loop_toggle.text = "Loop (5 s)"
	loop_toggle.toggled.connect(func(enabled: bool): loop_enabled_requested.emit(enabled))
	playback_controls.add_child(loop_toggle)
	timeline = HSlider.new()
	timeline.custom_minimum_size = Vector2(180.0, 20.0)
	timeline.size_flags_horizontal = Control.SIZE_EXPAND_FILL
	timeline.step = 1.0
	timeline.value_changed.connect(func(value: float): scrub_requested.emit(value))
	controls.add_child(timeline)

	view_panel = _panel()
	var camera_controls := _stack(view_panel, "VIEWPORTS")
	var camera_buttons := HFlowContainer.new()
	_flow_layout(camera_buttons)
	camera_controls.add_child(camera_buttons)
	for mode in [0, 1, 2]:
		var camera_button := _camera_button(mode)
		camera_buttons.add_child(camera_button)
		_view_buttons[mode] = camera_button
	_disabled_button(camera_buttons, "Pop out view", "A detachable camera viewport will be connected after multi-window rendering is implemented.")
	camera_controls.add_child(_label("Select one or more views. Drag with right mouse in Mouse Camera.", 10))

	monte_carlo_panel = _panel()
	var monte_carlo := _stack(monte_carlo_panel, "MONTE CARLO")
	monte_carlo.add_child(_label("Single-flight mode", 11))
	var monte_carlo_actions := HFlowContainer.new()
	_flow_layout(monte_carlo_actions)
	monte_carlo.add_child(monte_carlo_actions)
	_disabled_button(monte_carlo_actions, "◀ Flight", "Monte Carlo batch import is not available yet.")
	_disabled_button(monte_carlo_actions, "Focus", "Flight selection will be available after Monte Carlo import is implemented.")
	_disabled_button(monte_carlo_actions, "Flight ▶", "Monte Carlo batch import is not available yet.")

	graphs_tab = Button.new()
	graphs_tab.text = "G\nR\nA\nP\nH\nS"
	graphs_tab.tooltip_text = "Open graphs"
	graphs_tab.add_theme_stylebox_override("normal", _nav_style(Color("0d2b3e"), Color("3b6b82")))
	graphs_tab.add_theme_stylebox_override("hover", _nav_style(Color("174b6a"), Color("73d6ff")))
	graphs_tab.pressed.connect(_toggle_graph_drawer)
	add_child(graphs_tab)

	graphs_panel = _panel()
	graphs_panel.visible = false
	var analysis_stack := _stack(graphs_panel, "GRAPHS")
	var close_graphs := Button.new()
	close_graphs.text = "Close graphs"
	close_graphs.tooltip_text = "Close graph drawer"
	close_graphs.pressed.connect(_toggle_graph_drawer)
	analysis_stack.add_child(close_graphs)
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
		graph.title = TelemetrySchema.display_name_with_unit(column)
		graph.pop_out_requested.connect(_open_graph_window)
		graph.expanded_changed.connect(_move_graph_for_expansion)
		graph_list.add_child(graph)
		graph_widgets.append(graph)
	var collapsed := ScrollContainer.new()
	collapsed.name = "Collapsed"
	collapsed.size_flags_vertical = Control.SIZE_EXPAND_FILL
	tabs.add_child(collapsed)
	collapsed_graph_list = VBoxContainer.new()
	collapsed_graph_list.size_flags_horizontal = Control.SIZE_EXPAND_FILL
	collapsed.add_child(collapsed_graph_list)
	var body := VBoxContainer.new()
	body.name = "Body Frame"
	body.add_child(_label("CG & tank frame", 18))
	body.add_child(_label("Cyan: CG\nGold: CP\nBlue: fuel\nGreen: oxidizer\nOrange: thrust\nGreen: velocity", 14))
	tabs.add_child(body)

	file_dialog = FileDialog.new()
	file_dialog.file_mode = FileDialog.FILE_MODE_OPEN_FILE
	file_dialog.access = FileDialog.ACCESS_FILESYSTEM
	file_dialog.filters = PackedStringArray(["*.csv ; Telemetry CSV"])
	file_dialog.file_selected.connect(_file_selected)
	add_child(file_dialog)

func _create_viewports() -> void:
	if not is_instance_valid(flight_world) or _viewports_ready:
		return
	for mode in [0, 1, 2]:
		var pane := PanelContainer.new()
		pane.name = "%s Pane" % VIEW_NAMES[mode]
		pane.add_theme_stylebox_override("panel", _pane_style())
		view_grid.add_child(pane)
		var pane_content := Control.new()
		pane_content.set_anchors_and_offsets_preset(Control.PRESET_FULL_RECT)
		pane.add_child(pane_content)
		var container := SubViewportContainer.new()
		container.stretch = true
		container.mouse_filter = Control.MOUSE_FILTER_STOP
		container.set_anchors_and_offsets_preset(Control.PRESET_FULL_RECT)
		pane_content.add_child(container)
		var viewport := SubViewport.new()
		viewport.size = Vector2i(960, 540)
		viewport.render_target_update_mode = SubViewport.UPDATE_ALWAYS
		viewport.handle_input_locally = false
		viewport.world_3d = get_viewport().world_3d
		container.add_child(viewport)
		var camera := Camera3D.new()
		camera.fov = 62.0
		camera.current = true
		viewport.add_child(camera)
		var label := Label.new()
		label.text = VIEW_NAMES[mode]
		label.position = Vector2(9.0, 7.0)
		label.mouse_filter = Control.MOUSE_FILTER_IGNORE
		label.add_theme_font_size_override("font_size", 11)
		label.add_theme_color_override("font_color", Color("e9f7ff"))
		pane_content.add_child(label)
		container.gui_input.connect(func(event: InputEvent): _viewport_gui_input(mode, event))
		_view_cameras[mode] = camera
		_view_panes[mode] = pane
	_viewports_ready = true
	_layout_ui()

func _viewport_gui_input(mode: int, event: InputEvent) -> void:
	if mode == 1 and event is InputEventMouseMotion and Input.is_mouse_button_pressed(MOUSE_BUTTON_RIGHT):
		orbit_requested.emit((event as InputEventMouseMotion).relative)

func _choose_telemetry() -> void:
	_file_action = "telemetry"
	file_dialog.popup_centered_ratio(0.75)

func _choose_reference() -> void:
	_file_action = "reference"
	file_dialog.popup_centered_ratio(0.75)

func _request_terrain() -> void:
	var latitude_text := terrain_latitude_field.text.strip_edges()
	var longitude_text := terrain_longitude_field.text.strip_edges()
	if not latitude_text.is_valid_float() or not longitude_text.is_valid_float():
		show_status("Enter numeric latitude and longitude before loading terrain.")
		return
	var api_key := terrain_key_field.text.strip_edges()
	if api_key.is_empty():
		show_status("Enter an OpenTopography API key before loading terrain.")
		return
	terrain_load_requested.emit(latitude_text.to_float(), longitude_text.to_float(), api_key, float(terrain_area_selector.get_selected_id()))
	terrain_key_field.clear()

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
	window.transient = true
	window.unresizable = false
	get_window().gui_embed_subwindows = false
	get_tree().current_scene.add_child(window)
	var graph := TelemetryGraphScript.new()
	graph.database = database
	graph.column = source.column
	graph.title = source.title
	graph.line_color = source.line_color
	graph.show_expand_toggle = false
	graph.set_anchors_and_offsets_preset(Control.PRESET_FULL_RECT)
	window.add_child(graph)
	popup_graph_widgets.append(graph)
	window.close_requested.connect(func(): popup_graph_widgets.erase(graph); window.queue_free())
	window.popup_centered()

func _move_graph_for_expansion(graph: Control, expanded: bool) -> void:
	var target_list := graph_list if expanded else collapsed_graph_list
	graph.reparent(target_list)
	var target_index := 0
	for sibling in target_list.get_children():
		if sibling != graph and sibling.column < graph.column:
			target_index += 1
	target_list.move_child(graph, target_index)

func _toggle_dropdown(target: Control) -> void:
	var show_target := not target.visible
	force_panel.hide()
	reference_panel.hide()
	if show_target:
		target.show()
	_layout_ui()

func _toggle_graph_drawer() -> void:
	_graphs_open = not _graphs_open
	graphs_panel.visible = _graphs_open
	graphs_tab.visible = not _graphs_open
	_layout_ui()

func _toggle_camera(mode: int) -> void:
	if mode in _active_view_modes:
		if _active_view_modes.size() == 1:
			return
		_active_view_modes.erase(mode)
		if _primary_view_mode == mode:
			_primary_view_mode = _active_view_modes[0]
	else:
		_active_view_modes.append(mode)
		_primary_view_mode = mode
	_update_camera_button_styles()
	_layout_camera_panes()

func _update_camera_button_styles() -> void:
	for mode in _view_buttons:
		var button := _view_buttons[mode] as Button
		var selected: bool = mode in _active_view_modes
		button.add_theme_stylebox_override("normal", _camera_button_style(selected, false))
		button.add_theme_stylebox_override("hover", _camera_button_style(selected, true))
		button.tooltip_text = "Hide %s" % VIEW_NAMES[mode] if selected else "Show %s" % VIEW_NAMES[mode]

func _panel(translucent: bool = false) -> PanelContainer:
	var panel := PanelContainer.new()
	panel.add_theme_stylebox_override("panel", _panel_style(translucent))
	add_child(panel)
	return panel

func _stack(panel: PanelContainer, heading: String) -> VBoxContainer:
	var stack := VBoxContainer.new()
	stack.add_theme_constant_override("separation", 6)
	panel.add_child(stack)
	var handle := _label(heading, 11)
	handle.add_theme_color_override("font_color", Color("7fb5cf"))
	stack.add_child(handle)
	return stack

func _navigation_button(parent: Container, text: String, action: Callable = Callable(), tooltip: String = "") -> void:
	var button := Button.new()
	button.text = text
	button.custom_minimum_size = Vector2(92.0, 30.0)
	_apply_navigation_style(button)
	if action.is_valid():
		button.pressed.connect(action)
	else:
		button.disabled = true
		button.tooltip_text = tooltip
	parent.add_child(button)

func _camera_button(mode: int) -> Button:
	var button := Button.new()
	button.text = VIEW_NAMES[mode]
	button.custom_minimum_size = Vector2(96.0, 30.0)
	button.pressed.connect(func(): _toggle_camera(mode))
	return button

func _button(parent: Container, text: String, action: Callable) -> void:
	var button := Button.new()
	button.text = text
	button.pressed.connect(action)
	parent.add_child(button)

func _disabled_button(parent: Container, text: String, tooltip: String) -> void:
	var button := Button.new()
	button.text = text
	button.disabled = true
	button.tooltip_text = tooltip
	parent.add_child(button)

func _label(text: String, font_size: int) -> Label:
	var label := Label.new()
	label.text = text
	label.add_theme_font_size_override("font_size", font_size)
	label.autowrap_mode = TextServer.AUTOWRAP_WORD_SMART
	return label

func _terrain_field(placeholder: String, secret: bool) -> LineEdit:
	var field := LineEdit.new()
	field.placeholder_text = placeholder
	field.custom_minimum_size = Vector2(116.0 if not secret else 170.0, 30.0)
	field.secret = secret
	return field

func _flow_layout(container: FlowContainer) -> void:
	container.add_theme_constant_override("h_separation", 6)
	container.add_theme_constant_override("v_separation", 4)

func _layout_ui() -> void:
	if size.x <= 1.0 or size.y <= 1.0:
		return
	var margin := clampf(minf(size.x, size.y) * 0.012, 8.0, 16.0)
	var compact := size.y < 800.0
	var navigation_height := 46.0
	var dropdown_height := 172.0 if reference_panel.visible else (106.0 if force_panel.visible else 0.0)
	var top_height := navigation_height + dropdown_height
	var bottom_height := clampf(size.y * (0.245 if compact else 0.215), 152.0, 220.0)
	_layout_panel(top_bar, Vector2.ZERO, Vector2(size.x, navigation_height))
	top_dropdown_band.visible = dropdown_height > 0.0
	if top_dropdown_band.visible:
		_layout_panel(top_dropdown_band, Vector2(0.0, navigation_height), Vector2(size.x, dropdown_height))
	if force_panel.visible:
		_layout_panel(force_panel, Vector2(margin, navigation_height), Vector2(minf(400.0, size.x - margin * 2.0), dropdown_height))
	if reference_panel.visible:
		_layout_panel(reference_panel, Vector2(margin + 104.0, navigation_height), Vector2(minf(630.0, size.x - margin * 2.0), dropdown_height))
	var bottom_y := size.y - bottom_height
	_layout_panel(bottom_band, Vector2(0.0, bottom_y), Vector2(size.x, bottom_height))
	var playback_width := clampf(size.x * 0.31, 360.0, 500.0)
	var view_width := clampf(size.x * 0.23, 300.0, 400.0)
	var monte_width := clampf(size.x * 0.18, 235.0, 300.0)
	var controls_gap := 10.0
	var controls_width := playback_width + view_width + monte_width + controls_gap * 2.0
	var controls_left := maxf(0.0, (size.x - controls_width) * 0.5)
	_layout_panel(playback_panel, Vector2(controls_left, bottom_y), Vector2(playback_width, bottom_height))
	_layout_panel(view_panel, Vector2(controls_left + playback_width + controls_gap, bottom_y), Vector2(view_width, bottom_height))
	_layout_panel(monte_carlo_panel, Vector2(controls_left + playback_width + view_width + controls_gap * 2.0, bottom_y), Vector2(monte_width, bottom_height))

	var workspace_top := top_height + margin
	var workspace_height := maxf(80.0, bottom_y - workspace_top - margin)
	var drawer_width := clampf(size.x * 0.29, 315.0, 470.0)
	var tab_width := 34.0
	var workspace_right := size.x - margin - (drawer_width + margin if _graphs_open else tab_width)
	workspace.position = Vector2(margin, workspace_top)
	workspace.size = Vector2(maxf(140.0, workspace_right - margin), workspace_height)
	view_grid.size = workspace.size
	_layout_panel(telemetry_panel, workspace.position + Vector2(10.0, 10.0), Vector2(minf(248.0, workspace.size.x - 20.0), 88.0))
	if _graphs_open:
		_layout_panel(graphs_panel, Vector2(size.x - drawer_width, workspace_top), Vector2(drawer_width, workspace_height))
	else:
		graphs_tab.position = Vector2(size.x - tab_width, workspace_top)
		graphs_tab.size = Vector2(tab_width, workspace_height)
	graph_list.custom_minimum_size.x = maxf(250.0, graphs_panel.size.x - 30.0)
	collapsed_graph_list.custom_minimum_size.x = maxf(250.0, graphs_panel.size.x - 30.0)
	for graph in graph_widgets:
		graph.set_card_width(graphs_panel.size.x - 30.0)
	_update_camera_button_styles()
	_layout_camera_panes()

func _layout_camera_panes() -> void:
	if not _viewports_ready:
		return
	for mode in _view_panes:
		(_view_panes[mode] as Control).visible = mode in _active_view_modes
	var active := _active_view_modes.duplicate()
	active.erase(_primary_view_mode)
	active.push_front(_primary_view_mode)
	var gap := 5.0
	var bounds := Rect2(Vector2.ZERO, view_grid.size)
	if active.size() == 1:
		_place_pane(active[0], bounds)
	elif active.size() == 2:
		var half_width := (bounds.size.x - gap) * 0.5
		_place_pane(active[0], Rect2(bounds.position, Vector2(half_width, bounds.size.y)))
		_place_pane(active[1], Rect2(Vector2(half_width + gap, 0.0), Vector2(half_width, bounds.size.y)))
	else:
		var primary_width := (bounds.size.x - gap) * 0.62
		var secondary_width := bounds.size.x - primary_width - gap
		var secondary_height := (bounds.size.y - gap) * 0.5
		_place_pane(active[0], Rect2(Vector2.ZERO, Vector2(primary_width, bounds.size.y)))
		_place_pane(active[1], Rect2(Vector2(primary_width + gap, 0.0), Vector2(secondary_width, secondary_height)))
		_place_pane(active[2], Rect2(Vector2(primary_width + gap, secondary_height + gap), Vector2(secondary_width, secondary_height)))

func _place_pane(mode: int, rect: Rect2) -> void:
	var pane := _view_panes[mode] as Control
	pane.position = rect.position
	pane.size = rect.size

func _layout_panel(panel: Control, panel_position: Vector2, panel_size: Vector2) -> void:
	panel.position = panel_position
	panel.size = panel_size

func _panel_style(translucent: bool = false) -> StyleBoxFlat:
	var style := StyleBoxFlat.new()
	var base: Color = _presentation_colors.get("ui_base", Color("062e49"))
	var highlight: Color = _presentation_colors.get("ui_highlight", Color("66cbed"))
	style.bg_color = Color(base.darkened(0.54), 0.79) if translucent else Color(base.darkened(0.57), 0.96)
	style.border_color = highlight.darkened(0.48)
	style.set_border_width_all(1)
	style.corner_radius_top_left = 7
	style.corner_radius_top_right = 7
	style.corner_radius_bottom_left = 7
	style.corner_radius_bottom_right = 7
	style.content_margin_left = 10
	style.content_margin_right = 10
	style.content_margin_top = 7
	style.content_margin_bottom = 7
	return style

func _pane_style() -> StyleBoxFlat:
	var style := StyleBoxFlat.new()
	style.bg_color = Color("02070c")
	style.border_color = (_presentation_colors.get("ui_highlight", Color("66cbed")) as Color).darkened(0.32)
	style.set_border_width_all(1)
	style.corner_radius_top_left = 4
	style.corner_radius_top_right = 4
	style.corner_radius_bottom_left = 4
	style.corner_radius_bottom_right = 4
	return style

func _nav_style(fill: Color, border: Color) -> StyleBoxFlat:
	var style := StyleBoxFlat.new()
	style.bg_color = fill
	style.border_color = border
	style.set_border_width_all(1)
	style.corner_radius_top_left = 5
	style.corner_radius_top_right = 5
	style.corner_radius_bottom_left = 5
	style.corner_radius_bottom_right = 5
	return style

func _camera_button_style(selected: bool, hover: bool) -> StyleBoxFlat:
	var base: Color = _presentation_colors.get("ui_base", Color("062e49"))
	var highlight: Color = _presentation_colors.get("ui_highlight", Color("66cbed"))
	var fill := base.lightened(0.10) if selected else base.darkened(0.18)
	var border := highlight.lightened(0.12) if selected else highlight.darkened(0.48)
	if hover:
		fill = base.lightened(0.23) if selected else base.lightened(0.08)
	return _nav_style(fill, border)

func _apply_navigation_style(button: Button) -> void:
	var base: Color = _presentation_colors.get("ui_base", Color("062e49"))
	var highlight: Color = _presentation_colors.get("ui_highlight", Color("66cbed"))
	button.add_theme_stylebox_override("normal", _nav_style(base.darkened(0.18), highlight.darkened(0.48)))
	button.add_theme_stylebox_override("hover", _nav_style(base.lightened(0.12), highlight.lightened(0.08)))
	button.add_theme_stylebox_override("disabled", _nav_style(base.darkened(0.42), base.darkened(0.18)))
