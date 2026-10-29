## Full-screen, extensible application settings page.
extends Control

signal back_requested

var presentation: Node
var _backdrop: ColorRect
var _appearance_page: ScrollContainer
var _status: Label

func configure(presentation_settings: Node) -> void:
	presentation = presentation_settings

func _ready() -> void:
	set_anchors_and_offsets_preset(Control.PRESET_FULL_RECT)
	_build()
	if presentation != null:
		presentation.changed.connect(_apply_presentation)
		_apply_presentation(presentation.snapshot())

func _build() -> void:
	_backdrop = ColorRect.new()
	_backdrop.set_anchors_and_offsets_preset(Control.PRESET_FULL_RECT)
	add_child(_backdrop)
	var margin := MarginContainer.new()
	margin.set_anchors_and_offsets_preset(Control.PRESET_FULL_RECT)
	margin.add_theme_constant_override("margin_left", 28)
	margin.add_theme_constant_override("margin_top", 24)
	margin.add_theme_constant_override("margin_right", 28)
	margin.add_theme_constant_override("margin_bottom", 24)
	add_child(margin)
	var frame := PanelContainer.new()
	margin.add_child(frame)
	var content := VBoxContainer.new()
	content.add_theme_constant_override("separation", 14)
	frame.add_child(content)
	var header := HBoxContainer.new()
	content.add_child(header)
	var back := Button.new()
	back.text = "‹  MAIN MENU"
	back.custom_minimum_size = Vector2(144, 42)
	back.pressed.connect(func(): back_requested.emit())
	header.add_child(back)
	var title := Label.new()
	title.text = "SETTINGS"
	title.horizontal_alignment = HORIZONTAL_ALIGNMENT_CENTER
	title.add_theme_font_size_override("font_size", 28)
	title.size_flags_horizontal = Control.SIZE_EXPAND_FILL
	header.add_child(title)
	var spacer := Control.new()
	spacer.custom_minimum_size.x = 144
	header.add_child(spacer)
	var subtitle := Label.new()
	subtitle.text = "Customize the visualizer presentation.  Vehicle-specific values belong in Vehicle Selector."
	subtitle.horizontal_alignment = HORIZONTAL_ALIGNMENT_CENTER
	subtitle.autowrap_mode = TextServer.AUTOWRAP_WORD_SMART
	content.add_child(subtitle)
	var body := HBoxContainer.new()
	body.size_flags_vertical = Control.SIZE_EXPAND_FILL
	body.add_theme_constant_override("separation", 16)
	content.add_child(body)
	var categories := VBoxContainer.new()
	categories.custom_minimum_size.x = 190
	categories.add_theme_constant_override("separation", 8)
	body.add_child(categories)
	_add_category(categories, "APPEARANCE", true, "Colors for the UI and visual overlays.")
	_add_category(categories, "DISPLAY DEFAULTS", false, "Default displayed layers will be added after the v1 display plan is finalized.")
	_add_category(categories, "DEFAULT LOCATION", false, "Persistent terrain location will be added after the location workflow is finalized.")
	_add_category(categories, "CONSTANTS", false, "Non-vehicle numerical constants will be added with the technical specification.")
	_appearance_page = ScrollContainer.new()
	_appearance_page.size_flags_horizontal = Control.SIZE_EXPAND_FILL
	_appearance_page.size_flags_vertical = Control.SIZE_EXPAND_FILL
	body.add_child(_appearance_page)
	var appearance_panel := PanelContainer.new()
	appearance_panel.size_flags_horizontal = Control.SIZE_EXPAND_FILL
	_appearance_page.add_child(appearance_panel)
	var appearance := VBoxContainer.new()
	appearance.add_theme_constant_override("separation", 12)
	appearance_panel.add_child(appearance)
	var heading := Label.new()
	heading.text = "APPEARANCE"
	heading.add_theme_font_size_override("font_size", 20)
	appearance.add_child(heading)
	var help := Label.new()
	help.text = "Changes apply immediately and are stored on this computer.  Terrain, sky, and imported model materials remain unchanged."
	help.autowrap_mode = TextServer.AUTOWRAP_WORD_SMART
	appearance.add_child(help)
	_add_color_group(appearance, "Interface", [["UI base", "ui_base"], ["UI highlight", "ui_highlight"]])
	_add_color_group(appearance, "Flight geometry", [["Trajectory", "trace"], ["Reference trajectory", "reference"], ["Ground projection", "projection"], ["Landing aids", "landing"], ["Uncertainty", "uncertainty"]])
	_add_color_group(appearance, "Vehicle overlays", [["Center of gravity", "cg"], ["Center of pressure", "cp"], ["Thrust", "thrust"], ["Velocity", "velocity"], ["Propellant", "propellant"], ["RCS", "rcs"]])
	_status = Label.new()
	_status.text = "Changes are saved automatically."
	_status.add_theme_color_override("font_color", Color("a9c7d9"))
	appearance.add_child(_status)

func _add_category(parent: Container, text: String, enabled: bool, tooltip: String) -> void:
	var button := Button.new()
	button.text = text
	button.custom_minimum_size = Vector2(180, 42)
	button.disabled = not enabled
	button.tooltip_text = tooltip
	parent.add_child(button)

func _add_color_group(parent: Container, title: String, definitions: Array) -> void:
	var section := VBoxContainer.new()
	section.add_theme_constant_override("separation", 7)
	parent.add_child(section)
	var heading := Label.new()
	heading.text = title.to_upper()
	heading.add_theme_font_size_override("font_size", 14)
	section.add_child(heading)
	var grid := GridContainer.new()
	grid.columns = 2
	grid.add_theme_constant_override("h_separation", 12)
	grid.add_theme_constant_override("v_separation", 6)
	section.add_child(grid)
	for definition in definitions:
		var label := Label.new()
		label.text = definition[0]
		label.vertical_alignment = VERTICAL_ALIGNMENT_CENTER
		grid.add_child(label)
		var picker := ColorPickerButton.new()
		picker.custom_minimum_size = Vector2(132, 30)
		picker.color = presentation.colors[definition[1]] if presentation != null else Color.WHITE
		picker.color_changed.connect(func(color: Color):
			if presentation != null:
				presentation.set_color(definition[1], color)
		)
		grid.add_child(picker)

func _apply_presentation(colors: Dictionary) -> void:
	_backdrop.color = (colors["ui_base"] as Color).darkened(0.68)
	if presentation != null:
		presentation.apply_controls(self)
