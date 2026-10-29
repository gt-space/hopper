## Full-screen vehicle carousel and persistent profile editor.
extends Control

signal back_requested
signal vehicle_selected(profile: Dictionary)

class StarfieldBackdrop extends Control:
	func _ready() -> void:
		mouse_filter = Control.MOUSE_FILTER_IGNORE
		resized.connect(queue_redraw)
		queue_redraw()

	func _draw() -> void:
		draw_rect(Rect2(Vector2.ZERO, size), Color("020714"))
		for index in 180:
			var seed := float(index) * 1.618
			var position := Vector2(_random(seed) * size.x, _random(seed + 23.0) * size.y)
			var brightness := 0.28 + _random(seed + 47.0) * 0.62
			var radius := 0.55 + _random(seed + 71.0) * 1.35
			var tint := Color("a9dcff") if _random(seed + 97.0) > 0.16 else Color("c5a8ff")
			tint.a = brightness
			draw_circle(position, radius, tint)
			if radius > 1.55:
				var flare := tint
				flare.a *= 0.32
				draw_line(position - Vector2(radius * 2.4, 0.0), position + Vector2(radius * 2.4, 0.0), flare, 0.6)
				draw_line(position - Vector2(0.0, radius * 2.4), position + Vector2(0.0, radius * 2.4), flare, 0.6)

	func _random(seed: float) -> float:
		return fposmod(sin(seed * 12.9898 + 78.233) * 43758.5453, 1.0)

var vehicle_store: Node
var presentation: Node
var _backdrop: ColorRect
var _preview_viewport: SubViewport
var _preview_container: SubViewportContainer
var _preview_model: Node3D
var _preview_side_models: Array[Node3D] = []
var _profile_name: Label
var _carousel_caption: Label
var _selection_button: Button
var _save_button: Button
var _status: Label
var _file_dialog: FileDialog
var _fields: Dictionary = {}
var _carousel_index := 0
var _spin_boost := 0.0

func configure(store: Node, presentation_settings: Node) -> void:
	vehicle_store = store
	presentation = presentation_settings

func _ready() -> void:
	set_anchors_and_offsets_preset(Control.PRESET_FULL_RECT)
	_build()
	if vehicle_store != null:
		vehicle_store.profiles_changed.connect(_profiles_changed)
		vehicle_store.selected_profile_changed.connect(func(_profile: Dictionary): _profiles_changed())
		_focus_selected_profile()
	if presentation != null:
		presentation.changed.connect(_apply_presentation)
		_apply_presentation(presentation.snapshot())

func _process(delta: float) -> void:
	if is_instance_valid(_preview_model):
		_preview_model.rotate_y(delta * (0.55 + _spin_boost * 2.8))
	for side_model in _preview_side_models:
		if is_instance_valid(side_model):
			side_model.rotate_y(delta * 0.24)
	_spin_boost = maxf(_spin_boost - delta * 1.8, 0.0)

func _build() -> void:
	_backdrop = ColorRect.new()
	_backdrop.set_anchors_and_offsets_preset(Control.PRESET_FULL_RECT)
	add_child(_backdrop)
	var margin := MarginContainer.new()
	margin.set_anchors_and_offsets_preset(Control.PRESET_FULL_RECT)
	margin.add_theme_constant_override("margin_left", 24)
	margin.add_theme_constant_override("margin_top", 20)
	margin.add_theme_constant_override("margin_right", 24)
	margin.add_theme_constant_override("margin_bottom", 20)
	add_child(margin)
	var frame := PanelContainer.new()
	margin.add_child(frame)
	var stack := VBoxContainer.new()
	stack.add_theme_constant_override("separation", 12)
	frame.add_child(stack)
	var header := HBoxContainer.new()
	stack.add_child(header)
	var exit := Button.new()
	exit.text = "‹  EXIT"
	exit.custom_minimum_size = Vector2(118, 42)
	exit.pressed.connect(func(): back_requested.emit())
	header.add_child(exit)
	var title := Label.new()
	title.text = "VEHICLE SELECTOR"
	title.horizontal_alignment = HORIZONTAL_ALIGNMENT_CENTER
	title.add_theme_font_size_override("font_size", 27)
	title.size_flags_horizontal = Control.SIZE_EXPAND_FILL
	header.add_child(title)
	var spacer := Control.new()
	spacer.custom_minimum_size.x = 118
	header.add_child(spacer)
	# The selector follows the reference sketch: carousel first, actions below,
	# then the editable vehicle profile in its own lower region.
	var carousel_panel := PanelContainer.new()
	carousel_panel.custom_minimum_size = Vector2(0, 260)
	carousel_panel.size_flags_horizontal = Control.SIZE_EXPAND_FILL
	carousel_panel.size_flags_vertical = Control.SIZE_EXPAND_FILL
	stack.add_child(carousel_panel)
	var carousel := VBoxContainer.new()
	carousel.add_theme_constant_override("separation", 8)
	carousel_panel.add_child(carousel)
	var carousel_heading := Label.new()
	carousel_heading.text = "VEHICLE CAROUSEL"
	carousel_heading.horizontal_alignment = HORIZONTAL_ALIGNMENT_CENTER
	carousel_heading.add_theme_font_size_override("font_size", 15)
	carousel.add_child(carousel_heading)
	var preview_stack := Control.new()
	preview_stack.custom_minimum_size = Vector2(640, 180)
	preview_stack.size_flags_vertical = Control.SIZE_EXPAND_FILL
	preview_stack.clip_contents = true
	carousel.add_child(preview_stack)
	var starfield := StarfieldBackdrop.new()
	starfield.set_anchors_and_offsets_preset(Control.PRESET_FULL_RECT)
	preview_stack.add_child(starfield)
	_preview_container = SubViewportContainer.new()
	_preview_container.stretch = true
	_preview_container.set_anchors_and_offsets_preset(Control.PRESET_FULL_RECT)
	preview_stack.add_child(_preview_container)
	_profile_name = Label.new()
	_profile_name.horizontal_alignment = HORIZONTAL_ALIGNMENT_CENTER
	_profile_name.add_theme_font_size_override("font_size", 20)
	carousel.add_child(_profile_name)
	_carousel_caption = Label.new()
	_carousel_caption.horizontal_alignment = HORIZONTAL_ALIGNMENT_CENTER
	_carousel_caption.autowrap_mode = TextServer.AUTOWRAP_WORD_SMART
	_carousel_caption.add_theme_font_size_override("font_size", 11)
	carousel.add_child(_carousel_caption)
	var navigation := HBoxContainer.new()
	navigation.alignment = BoxContainer.ALIGNMENT_CENTER
	navigation.add_theme_constant_override("separation", 12)
	carousel.add_child(navigation)
	var previous := Button.new()
	previous.text = "◀  PREVIOUS"
	previous.custom_minimum_size = Vector2(140, 42)
	previous.pressed.connect(func(): _advance_carousel(-1))
	navigation.add_child(previous)
	var next := Button.new()
	next.text = "NEXT  ▶"
	next.custom_minimum_size = Vector2(140, 42)
	next.pressed.connect(func(): _advance_carousel(1))
	navigation.add_child(next)
	var actions := HBoxContainer.new()
	actions.alignment = BoxContainer.ALIGNMENT_CENTER
	actions.add_theme_constant_override("separation", 8)
	stack.add_child(actions)
	var create := Button.new()
	create.text = "CREATE"
	create.custom_minimum_size = Vector2(160, 48)
	create.pressed.connect(_choose_glb)
	actions.add_child(create)
	_selection_button = Button.new()
	_selection_button.custom_minimum_size = Vector2(160, 48)
	_selection_button.pressed.connect(_select_current)
	actions.add_child(_selection_button)
	_save_button = Button.new()
	_save_button.text = "SAVE"
	_save_button.custom_minimum_size = Vector2(160, 48)
	_save_button.pressed.connect(_save_current)
	actions.add_child(_save_button)
	var form_scroll := ScrollContainer.new()
	form_scroll.custom_minimum_size = Vector2(0, 190)
	form_scroll.size_flags_horizontal = Control.SIZE_EXPAND_FILL
	form_scroll.size_flags_vertical = Control.SIZE_EXPAND_FILL
	stack.add_child(form_scroll)
	var form_panel := PanelContainer.new()
	form_panel.size_flags_horizontal = Control.SIZE_EXPAND_FILL
	form_scroll.add_child(form_panel)
	var form := VBoxContainer.new()
	form.add_theme_constant_override("separation", 9)
	form_panel.add_child(form)
	var form_heading := Label.new()
	form_heading.text = "VEHICLE PROFILE"
	form_heading.add_theme_font_size_override("font_size", 17)
	form.add_child(form_heading)
	_add_text_field(form, "Name", "display_name")
	_add_number_group(form, "Dimensions (m)", [["Height", "height"], ["Width", "width"], ["Length", "length"]], "dimensions")
	_add_number_group(form, "RCS offset (m)", [["X", "x"], ["Y", "y"], ["Z", "z"]], "rcs_offset")
	_add_number_group(form, "TVC offset (m)", [["X", "x"], ["Y", "y"], ["Z", "z"]], "tvc_offset")
	_add_number_group(form, "Model transform", [["Scale", "model_scale"], ["Rotation X°", "rotation_x"], ["Rotation Y°", "rotation_y"], ["Rotation Z°", "rotation_z"]], "transform")
	_status = Label.new()
	_status.autowrap_mode = TextServer.AUTOWRAP_WORD_SMART
	_status.add_theme_color_override("font_color", Color("a9c7d9"))
	form.add_child(_status)
	_file_dialog = FileDialog.new()
	_file_dialog.access = FileDialog.ACCESS_FILESYSTEM
	_file_dialog.file_mode = FileDialog.FILE_MODE_OPEN_FILE
	_file_dialog.filters = PackedStringArray(["*.glb ; Binary glTF vehicle model"])
	_file_dialog.file_selected.connect(_create_from_glb)
	add_child(_file_dialog)
	_build_preview()

func _build_preview() -> void:
	_preview_viewport = SubViewport.new()
	_preview_viewport.own_world_3d = true
	_preview_viewport.transparent_bg = true
	_preview_viewport.render_target_update_mode = SubViewport.UPDATE_ALWAYS
	_preview_container.add_child(_preview_viewport)
	var environment := WorldEnvironment.new()
	var settings := Environment.new()
	settings.background_mode = Environment.BG_COLOR
	settings.background_color = Color(0.0, 0.0, 0.0, 0.0)
	settings.ambient_light_source = Environment.AMBIENT_SOURCE_COLOR
	settings.ambient_light_color = Color("a2c8db")
	settings.ambient_light_energy = 0.62
	environment.environment = settings
	_preview_viewport.add_child(environment)
	var light := DirectionalLight3D.new()
	light.rotation_degrees = Vector3(-42.0, -28.0, 0.0)
	light.light_energy = 2.0
	_preview_viewport.add_child(light)
	var camera := Camera3D.new()
	camera.look_at_from_position(Vector3(0.0, 3.0, 10.5), Vector3(0.0, 0.65, 0.0))
	camera.current = true
	_preview_viewport.add_child(camera)

func _add_text_field(parent: Container, label_text: String, key: String) -> void:
	var label := Label.new()
	label.text = label_text.to_upper()
	parent.add_child(label)
	var field := LineEdit.new()
	field.custom_minimum_size.y = 32
	parent.add_child(field)
	_fields[key] = field

func _add_number_group(parent: Container, title: String, definitions: Array, group: String) -> void:
	var label := Label.new()
	label.text = title.to_upper()
	parent.add_child(label)
	var grid := GridContainer.new()
	grid.columns = 2
	grid.add_theme_constant_override("h_separation", 8)
	grid.add_theme_constant_override("v_separation", 5)
	parent.add_child(grid)
	for definition in definitions:
		var name := Label.new()
		name.text = definition[0]
		name.vertical_alignment = VERTICAL_ALIGNMENT_CENTER
		grid.add_child(name)
		var field := SpinBox.new()
		field.min_value = 0.001 if definition[1] in ["height", "width", "length", "model_scale"] else -1000.0
		field.max_value = 1000.0
		field.step = 0.01
		field.custom_minimum_size = Vector2(122, 30)
		grid.add_child(field)
		_fields["%s.%s" % [group, definition[1]]] = field

func _profiles_changed() -> void:
	_focus_selected_profile()

func _focus_selected_profile() -> void:
	if vehicle_store == null or vehicle_store.profiles.is_empty():
		return
	for index in vehicle_store.profiles.size():
		if vehicle_store.profiles[index].get("id", "") == vehicle_store.selected_id:
			_carousel_index = index
			break
	_show_current_profile()

func _advance_carousel(direction: int) -> void:
	if vehicle_store == null or vehicle_store.profiles.is_empty():
		return
	_carousel_index = posmod(_carousel_index + direction, vehicle_store.profiles.size())
	_show_current_profile()

func _show_current_profile() -> void:
	if vehicle_store == null or vehicle_store.profiles.is_empty():
		return
	var profile: Dictionary = vehicle_store.profiles[_carousel_index]
	_profile_name.text = str(profile["display_name"])
	(_fields["display_name"] as LineEdit).text = str(profile["display_name"])
	_set_group_values("dimensions", profile.get("dimensions", {}))
	_set_group_values("rcs_offset", profile.get("rcs_offset", {}))
	_set_group_values("tvc_offset", profile.get("tvc_offset", {}))
	var transform: Dictionary = profile.get("model_rotation", {})
	(_fields["transform.model_scale"] as SpinBox).value = float(profile.get("model_scale", 1.0))
	(_fields["transform.rotation_x"] as SpinBox).value = float(transform.get("x", 0.0))
	(_fields["transform.rotation_y"] as SpinBox).value = float(transform.get("y", 0.0))
	(_fields["transform.rotation_z"] as SpinBox).value = float(transform.get("z", 0.0))
	_selection_button.text = "SELECTED" if profile.get("id", "") == vehicle_store.selected_id else "SELECT"
	_save_button.disabled = profile.get("id", "") == "hopper-default"
	_save_button.tooltip_text = "The built-in profile is read-only." if _save_button.disabled else ""
	_replace_preview(profile)
	_update_carousel_caption()
	_status.text = "Created %s" % Time.get_datetime_string_from_unix_time(int(profile.get("created_unix", 0))) if int(profile.get("created_unix", 0)) > 0 else "Built-in fallback vehicle"

func _set_group_values(group: String, values: Dictionary) -> void:
	for key in values:
		var field := _fields.get("%s.%s" % [group, key]) as SpinBox
		if field != null:
			field.value = float(values[key])

func _replace_preview(profile: Dictionary) -> void:
	if is_instance_valid(_preview_model):
		_preview_viewport.remove_child(_preview_model)
		_preview_model.queue_free()
	for side_model in _preview_side_models:
		if is_instance_valid(side_model):
			_preview_viewport.remove_child(side_model)
			side_model.queue_free()
	_preview_side_models.clear()
	_preview_model = vehicle_store.instantiate_model(profile) if vehicle_store != null else null
	if _preview_model != null:
		_preview_model.position = Vector3(0.0, 0.25, 0.0)
		_preview_model.scale *= 2.60
		_preview_viewport.add_child(_preview_model)
	if vehicle_store == null or vehicle_store.profiles.size() < 2:
		return
	var profile_count: int = vehicle_store.profiles.size()
	var previous_index := posmod(_carousel_index - 1, profile_count)
	var next_index := posmod(_carousel_index + 1, profile_count)
	if previous_index != _carousel_index:
		_add_side_preview(vehicle_store.profiles[previous_index], Vector3(-6.8, 0.15, 0.45), -20.0)
	if next_index != _carousel_index:
		_add_side_preview(vehicle_store.profiles[next_index], Vector3(6.8, 0.15, 0.45), 20.0)

func _add_side_preview(profile: Dictionary, position: Vector3, yaw: float) -> void:
	var side_model: Node3D = vehicle_store.instantiate_model(profile)
	if side_model == null:
		return
	side_model.position = position
	side_model.scale *= 1.36
	side_model.rotation_degrees += Vector3(0.0, yaw, 0.0)
	_preview_viewport.add_child(side_model)
	_preview_side_models.append(side_model)

func _update_carousel_caption() -> void:
	if vehicle_store == null or vehicle_store.profiles.size() < 2:
		_carousel_caption.text = "Add a GLB to begin building a vehicle carousel."
		return
	var count: int = vehicle_store.profiles.size()
	var previous_name := str(vehicle_store.profiles[posmod(_carousel_index - 1, count)].get("display_name", "Previous"))
	var next_name := str(vehicle_store.profiles[posmod(_carousel_index + 1, count)].get("display_name", "Next"))
	_carousel_caption.text = "‹  %s       %s  ›" % [previous_name, next_name]

func _choose_glb() -> void:
	_file_dialog.popup_centered_ratio(0.74)

func _create_from_glb(path: String) -> void:
	var profile: Dictionary = vehicle_store.create_profile_from_glb(path) if vehicle_store != null else {}
	if profile.is_empty():
		_status.text = "Could not create a profile from that GLB.  Keep the previous vehicle selected."
		return
	_status.text = "Created %s. Review the transform and press Save if needed." % profile["display_name"]
	_focus_selected_profile()

func _select_current() -> void:
	if vehicle_store == null or vehicle_store.profiles.is_empty():
		return
	var profile: Dictionary = vehicle_store.profiles[_carousel_index]
	if vehicle_store.select_profile(str(profile["id"])):
		_spin_boost = 1.0
		vehicle_selected.emit(vehicle_store.get_selected_profile())
		_show_current_profile()

func _save_current() -> void:
	if vehicle_store == null or vehicle_store.profiles.is_empty():
		return
	var profile: Dictionary = vehicle_store.profiles[_carousel_index].duplicate(true)
	if profile.get("id", "") == "hopper-default":
		return
	profile["display_name"] = (_fields["display_name"] as LineEdit).text.strip_edges()
	if str(profile["display_name"]).is_empty():
		_status.text = "A vehicle name is required."
		return
	profile["dimensions"] = _group_values("dimensions", ["height", "width", "length"])
	profile["rcs_offset"] = _group_values("rcs_offset", ["x", "y", "z"])
	profile["tvc_offset"] = _group_values("tvc_offset", ["x", "y", "z"])
	profile["model_scale"] = (_fields["transform.model_scale"] as SpinBox).value
	profile["model_rotation"] = {"x": (_fields["transform.rotation_x"] as SpinBox).value, "y": (_fields["transform.rotation_y"] as SpinBox).value, "z": (_fields["transform.rotation_z"] as SpinBox).value}
	if vehicle_store.save_profile(profile):
		_status.text = "Saved %s." % profile["display_name"]
		_replace_preview(profile)
		if profile.get("id", "") == vehicle_store.selected_id:
			vehicle_selected.emit(profile)
	else:
		_status.text = "Could not save this profile; the previous saved version is unchanged."

func _group_values(group: String, keys: Array) -> Dictionary:
	var values := {}
	for key in keys:
		values[key] = (_fields["%s.%s" % [group, key]] as SpinBox).value
	return values

func _apply_presentation(colors: Dictionary) -> void:
	_backdrop.color = (colors["ui_base"] as Color).darkened(0.68)
	if presentation != null:
		presentation.apply_controls(self)
