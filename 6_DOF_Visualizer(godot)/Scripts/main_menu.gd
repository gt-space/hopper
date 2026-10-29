## Front page for the visualizer. Profile selection and live-session setup are intentionally visual-only for now.
extends Control

signal load_mission_requested
signal settings_requested
signal vehicle_selector_requested

const ROCKET_SCENE = preload("res://Prefabs/Rocket.tscn")

var _preview_viewport: SubViewport
var _preview_container: SubViewportContainer
var _preview_rocket: Node3D
var _vehicle_card: PanelContainer
var _menu_stack: VBoxContainer
var _star_phase := 0.0
var _vehicle_name: Label
var vehicle_store: Node
var presentation: Node

func configure(store: Node, presentation_settings: Node) -> void:
	vehicle_store = store
	presentation = presentation_settings
	if vehicle_store != null:
		vehicle_store.selected_profile_changed.connect(_refresh_selected_vehicle)

func _ready() -> void:
	set_anchors_and_offsets_preset(Control.PRESET_FULL_RECT)
	_build()
	resized.connect(_layout)
	call_deferred("_layout")

func _process(delta: float) -> void:
	_star_phase += delta * 0.08
	if is_instance_valid(_preview_rocket):
		_preview_rocket.rotate_y(delta * 0.55)
	queue_redraw()

func _draw() -> void:
	draw_rect(Rect2(Vector2.ZERO, size), Color("040a13"))
	for star_index in range(78):
		var angle := float(star_index) * 2.39996 + _star_phase
		var radius := 25.0 + fmod(float(star_index) * 41.0, maxf(size.length() * 0.72, 1.0))
		var center := size * 0.5 + Vector2(cos(angle), sin(angle) * 0.62) * radius
		var brightness := 0.28 + fmod(float(star_index) * 0.37, 0.48)
		draw_circle(center, 0.7 + fmod(float(star_index), 2.0) * 0.35, Color(0.60, 0.82, 1.0, brightness))
	var glow := Rect2(size * 0.5 - Vector2(300.0, 210.0), Vector2(600.0, 420.0))
	draw_style_box(_glow_style(), glow)

func _build() -> void:
	var heading := Label.new()
	heading.text = "HOPPER"
	heading.horizontal_alignment = HORIZONTAL_ALIGNMENT_CENTER
	heading.add_theme_font_size_override("font_size", 46)
	heading.add_theme_color_override("font_color", Color("d8f4ff"))
	add_child(heading)
	var subtitle := Label.new()
	subtitle.text = "6-DOF TELEMETRY VISUALIZER"
	subtitle.horizontal_alignment = HORIZONTAL_ALIGNMENT_CENTER
	subtitle.add_theme_font_size_override("font_size", 16)
	subtitle.add_theme_color_override("font_color", Color("75bedf"))
	add_child(subtitle)
	_menu_stack = VBoxContainer.new()
	_menu_stack.alignment = BoxContainer.ALIGNMENT_CENTER
	_menu_stack.add_theme_constant_override("separation", 16)
	add_child(_menu_stack)
	_menu_stack.add_child(_menu_button("QUIT", func(): get_tree().quit()))
	_menu_stack.add_child(_menu_button("SETTINGS", func(): settings_requested.emit()))
	var mission_row := HBoxContainer.new()
	mission_row.alignment = BoxContainer.ALIGNMENT_CENTER
	mission_row.add_theme_constant_override("separation", 18)
	mission_row.add_child(_menu_button("LOAD MISSION", func(): load_mission_requested.emit()))
	mission_row.add_child(_menu_button("LIVE SESSION", Callable(), "Live-session setup will be connected after the telemetry protocol is finalized."))
	_menu_stack.add_child(mission_row)
	_vehicle_card = PanelContainer.new()
	_vehicle_card.add_theme_stylebox_override("panel", _panel_style())
	add_child(_vehicle_card)
	var vehicle_stack := VBoxContainer.new()
	vehicle_stack.add_theme_constant_override("separation", 6)
	_vehicle_card.add_child(vehicle_stack)
	var label := Label.new()
	label.text = "SELECTED VEHICLE"
	label.horizontal_alignment = HORIZONTAL_ALIGNMENT_CENTER
	label.add_theme_font_size_override("font_size", 12)
	label.add_theme_color_override("font_color", Color("87cce8"))
	vehicle_stack.add_child(label)
	_preview_container = SubViewportContainer.new()
	_preview_container.stretch = true
	_preview_container.custom_minimum_size = Vector2(276, 172)
	vehicle_stack.add_child(_preview_container)
	_vehicle_name = Label.new()
	_vehicle_name.text = "Hopper — Default Vehicle"
	_vehicle_name.horizontal_alignment = HORIZONTAL_ALIGNMENT_CENTER
	_vehicle_name.add_theme_font_size_override("font_size", 15)
	vehicle_stack.add_child(_vehicle_name)
	vehicle_stack.add_child(_menu_button("CHOOSE VEHICLE", func(): vehicle_selector_requested.emit()))
	_build_vehicle_preview()
	if presentation != null:
		presentation.changed.connect(apply_presentation)
		apply_presentation(presentation.snapshot())

func _build_vehicle_preview() -> void:
	_preview_viewport = SubViewport.new()
	# The preview must not inherit the flight world's 3D scene; otherwise it becomes a second rocket at the origin.
	_preview_viewport.own_world_3d = true
	_preview_viewport.transparent_bg = false
	_preview_viewport.render_target_update_mode = SubViewport.UPDATE_ALWAYS
	_preview_container.add_child(_preview_viewport)
	var environment := WorldEnvironment.new()
	var settings := Environment.new()
	settings.background_mode = Environment.BG_COLOR
	settings.background_color = Color("06111e")
	settings.ambient_light_source = Environment.AMBIENT_SOURCE_COLOR
	settings.ambient_light_color = Color("8cb9d0")
	settings.ambient_light_energy = 0.55
	environment.environment = settings
	_preview_viewport.add_child(environment)
	var light := DirectionalLight3D.new()
	light.rotation_degrees = Vector3(-38.0, -25.0, 0.0)
	light.light_energy = 2.1
	_preview_viewport.add_child(light)
	var profile: Dictionary = vehicle_store.get_selected_profile() if vehicle_store != null else {}
	_preview_rocket = vehicle_store.instantiate_model(profile) if not profile.is_empty() else ROCKET_SCENE.instantiate()
	_preview_rocket.rotation_degrees += Vector3(0.0, -24.0, 0.0)
	_preview_rocket.position.y = 0.18
	_preview_viewport.add_child(_preview_rocket)
	var camera := Camera3D.new()
	camera.look_at_from_position(Vector3(3.25, 1.85, 4.8), Vector3(0.0, 0.35, 0.0))
	camera.current = true
	_preview_viewport.add_child(camera)

func _refresh_selected_vehicle(profile: Dictionary) -> void:
	if is_instance_valid(_vehicle_name):
		_vehicle_name.text = str(profile.get("display_name", "Hopper — Default Vehicle"))
	if is_instance_valid(_preview_rocket):
		_preview_viewport.remove_child(_preview_rocket)
		_preview_rocket.queue_free()
		_preview_rocket = vehicle_store.instantiate_model(profile)
		if _preview_rocket != null:
			_preview_rocket.rotation_degrees += Vector3(0.0, -24.0, 0.0)
			_preview_rocket.position.y = 0.18
			_preview_viewport.add_child(_preview_rocket)

func apply_presentation(colors: Dictionary) -> void:
	if presentation == null:
		return
	presentation.apply_controls(self)
	if is_instance_valid(_vehicle_card):
		_vehicle_card.add_theme_stylebox_override("panel", _panel_style())

func _layout() -> void:
	if size.x < 2.0 or size.y < 2.0:
		return
	var compact := size.x < 900.0 or size.y < 650.0
	var menu_width := 472.0 if not compact else minf(414.0, size.x - 36.0)
	var menu_height := 252.0
	_menu_stack.position = Vector2((size.x - menu_width) * 0.5, (size.y - menu_height) * 0.5 - 4.0)
	_menu_stack.size = Vector2(menu_width, menu_height)
	get_child(0).position = Vector2(0.0, maxf(28.0, size.y * 0.13))
	get_child(0).size = Vector2(size.x, 58.0)
	get_child(1).position = Vector2(0.0, get_child(0).position.y + 54.0)
	get_child(1).size = Vector2(size.x, 28.0)
	var card_size := Vector2(345.0, 356.0) if not compact else Vector2(288.0, 322.0)
	_vehicle_card.size = card_size
	_vehicle_card.position = Vector2(size.x - card_size.x - 30.0, size.y - card_size.y - 28.0)

func _menu_button(text: String, action: Callable = Callable(), tooltip: String = "") -> Button:
	var button := Button.new()
	button.text = text
	button.custom_minimum_size = Vector2(205.0, 51.0)
	button.add_theme_stylebox_override("normal", _button_style(Color("66cbed"), Color("062e49")))
	button.add_theme_stylebox_override("hover", _button_style(Color("91e2ff"), Color("0b4569")))
	button.add_theme_stylebox_override("disabled", _button_style(Color("24465a"), Color("172633")))
	button.add_theme_color_override("font_color", Color("041420"))
	button.add_theme_color_override("font_hover_color", Color.WHITE)
	button.add_theme_color_override("font_disabled_color", Color("7290a0"))
	if action.is_valid():
		button.pressed.connect(action)
	else:
		button.disabled = true
		button.tooltip_text = tooltip
	return button

func _panel_style() -> StyleBoxFlat:
	var style := StyleBoxFlat.new()
	var base := Color("0b1c2b")
	var highlight := Color("4caed2")
	if presentation != null:
		base = presentation.colors["ui_base"]
		highlight = presentation.colors["ui_highlight"]
	style.bg_color = Color(base.darkened(0.30), 0.91)
	style.border_color = highlight
	style.set_border_width_all(1)
	style.corner_radius_top_left = 10
	style.corner_radius_top_right = 10
	style.corner_radius_bottom_left = 10
	style.corner_radius_bottom_right = 10
	style.content_margin_left = 12
	style.content_margin_right = 12
	style.content_margin_top = 10
	style.content_margin_bottom = 10
	return style

func _button_style(fill: Color, border: Color) -> StyleBoxFlat:
	var style := StyleBoxFlat.new()
	style.bg_color = fill
	style.border_color = border
	style.set_border_width_all(1)
	style.corner_radius_top_left = 6
	style.corner_radius_top_right = 6
	style.corner_radius_bottom_left = 6
	style.corner_radius_bottom_right = 6
	return style

func _glow_style() -> StyleBoxFlat:
	var style := StyleBoxFlat.new()
	style.bg_color = Color("1a4f7030")
	style.corner_radius_top_left = 260
	style.corner_radius_top_right = 260
	style.corner_radius_bottom_left = 260
	style.corner_radius_bottom_right = 260
	return style
