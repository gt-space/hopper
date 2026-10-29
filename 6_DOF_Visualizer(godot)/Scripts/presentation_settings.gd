## Stores the user-editable visual palette.  It is owned by Main rather than autoloaded.
extends Node

signal changed(colors: Dictionary)

const SETTINGS_PATH := "user://settings.cfg"

const DEFAULT_COLORS := {
	"ui_base": Color("062e49"),
	"ui_highlight": Color("66cbed"),
	"trace": Color("58d6ff"),
	"reference": Color("ffcf5c"),
	"projection": Color("ff334f"),
	"landing": Color("f2b84b"),
	"uncertainty": Color("bd8cff"),
	"cg": Color("56e0ff"),
	"cp": Color("ffcf5c"),
	"thrust": Color("ff7b38"),
	"velocity": Color("57e389"),
	"propellant": Color("36a8ff"),
	"rcs": Color("66d9ff"),
}

var colors: Dictionary = DEFAULT_COLORS.duplicate(true)

func load_from_disk() -> void:
	colors = DEFAULT_COLORS.duplicate(true)
	var config := ConfigFile.new()
	if config.load(SETTINGS_PATH) == OK:
		for key in DEFAULT_COLORS:
			var value: Variant = config.get_value("appearance", key, DEFAULT_COLORS[key])
			if value is Color:
				colors[key] = value
	changed.emit(snapshot())

func set_color(key: String, color: Color) -> void:
	if not DEFAULT_COLORS.has(key):
		return
	colors[key] = color
	_save()
	changed.emit(snapshot())

func snapshot() -> Dictionary:
	return colors.duplicate(true)

func apply_controls(root: Node) -> void:
	var base: Color = colors["ui_base"]
	var highlight: Color = colors["ui_highlight"]
	for child in root.get_children():
		if child is Button:
			var button := child as Button
			button.add_theme_stylebox_override("normal", _button_style(highlight, base))
			button.add_theme_stylebox_override("hover", _button_style(base.lightened(0.10), highlight))
			button.add_theme_stylebox_override("disabled", _button_style(base.darkened(0.16), base.darkened(0.24)))
			button.add_theme_color_override("font_color", Color("041420"))
			button.add_theme_color_override("font_hover_color", Color.WHITE)
			button.add_theme_color_override("font_disabled_color", Color("7290a0"))
		elif child is PanelContainer:
			var panel := child as PanelContainer
			panel.add_theme_stylebox_override("panel", _panel_style(base, highlight))
		elif child is LineEdit or child is OptionButton or child is SpinBox:
			child.add_theme_color_override("font_color", Color("e9f7ff"))
			child.add_theme_color_override("font_placeholder_color", Color("9bb8c8"))
		apply_controls(child)

func _save() -> void:
	var config := ConfigFile.new()
	for key in colors:
		config.set_value("appearance", key, colors[key])
	config.save(SETTINGS_PATH)

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

func _panel_style(base: Color, highlight: Color) -> StyleBoxFlat:
	var style := StyleBoxFlat.new()
	style.bg_color = Color(base.darkened(0.46), 0.93)
	style.border_color = highlight.darkened(0.20)
	style.set_border_width_all(1)
	style.corner_radius_top_left = 9
	style.corner_radius_top_right = 9
	style.corner_radius_bottom_left = 9
	style.corner_radius_bottom_right = 9
	style.content_margin_left = 12
	style.content_margin_right = 12
	style.content_margin_top = 10
	style.content_margin_bottom = 10
	return style
