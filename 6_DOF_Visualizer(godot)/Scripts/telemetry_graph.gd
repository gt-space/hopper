## Reusable graph card with hover inspection, collapse control, and pop-out support.
extends Control

const TelemetrySchema = preload("res://Scripts/telemetry_schema.gd")

signal pop_out_requested(graph: Control)

var database: Node
var column := TelemetrySchema.Column.UP
var cursor := 0
var title := "Altitude"
var line_color := Color("58d6ff")
var hover_index := -1
var _expanded := true
var _expand_toggle: CheckButton
var _open_button: Button

func _ready() -> void:
	custom_minimum_size = Vector2(280, 150)
	mouse_filter = Control.MOUSE_FILTER_STOP
	_expand_toggle = CheckButton.new()
	_expand_toggle.text = "Expanded"
	_expand_toggle.button_pressed = true
	_expand_toggle.position = Vector2(185, 3)
	_expand_toggle.toggled.connect(_set_expanded)
	add_child(_expand_toggle)
	_open_button = Button.new()
	_open_button.text = "Open"
	_open_button.position = Vector2(292, 3)
	_open_button.pressed.connect(func(): pop_out_requested.emit(self))
	add_child(_open_button)
	mouse_exited.connect(_clear_hover)
	if database != null:
		database.data_changed.connect(queue_redraw)
	resized.connect(_place_header_controls)
	call_deferred("_place_header_controls")

func set_card_width(width: float) -> void:
	custom_minimum_size.x = maxf(width, 280.0)

func _place_header_controls() -> void:
	if _open_button == null or _expand_toggle == null:
		return
	_open_button.position.x = maxf(size.x - _open_button.size.x - 8.0, 215.0)
	_expand_toggle.position.x = maxf(_open_button.position.x - _expand_toggle.size.x - 8.0, 115.0)

func _set_expanded(expanded: bool) -> void:
	_expanded = expanded
	custom_minimum_size.y = 150.0 if expanded else 32.0
	_expand_toggle.text = "Expanded" if expanded else "Collapsed"
	queue_redraw()

func _gui_input(event: InputEvent) -> void:
	if not _expanded or not event is InputEventMouseMotion:
		return
	var plot := _plot_rect()
	var motion := event as InputEventMouseMotion
	var mouse: Vector2 = motion.position
	if plot.has_point(mouse) and database != null and database.rows.size() > 1:
		hover_index = clampi(roundi((mouse.x - plot.position.x) / plot.size.x * float(database.rows.size() - 1)), 0, database.rows.size() - 1)
	else:
		hover_index = -1
	queue_redraw()

func _clear_hover() -> void:
	hover_index = -1
	queue_redraw()

func _draw() -> void:
	var rect := Rect2(Vector2.ZERO, size).grow(-3.0)
	draw_style_box(_panel_style(), rect)
	draw_string(ThemeDB.fallback_font, Vector2(9, 20), title, HORIZONTAL_ALIGNMENT_LEFT, maxf(_expand_toggle.position.x - 16.0, 100.0), 14, Color.WHITE)
	if not _expanded:
		return
	if database == null or database.rows.size() < 2:
		draw_string(ThemeDB.fallback_font, rect.get_center() + Vector2(-46, 5), "Waiting for data", HORIZONTAL_ALIGNMENT_LEFT, -1, 12, Color("a9b9c6"))
		return
	var plot := _plot_rect()
	var limits: Vector2 = database.value_range(column)
	var span := maxf(limits.y - limits.x, 0.0001)
	for i in range(1, 5):
		var y := plot.position.y + plot.size.y * float(i) / 5.0
		draw_line(Vector2(plot.position.x, y), Vector2(plot.end.x, y), Color("223848"), 1.0)
	var points := PackedVector2Array()
	var stride := maxi(1, ceili(float(database.rows.size()) / maxf(plot.size.x, 1.0)))
	for i in range(0, database.rows.size(), stride):
		points.append(_point_for(i, plot, limits, span))
	if points.size() > 1:
		draw_polyline(points, line_color, 1.8, true)
	var current_x := _point_for(cursor, plot, limits, span).x
	draw_line(Vector2(current_x, plot.position.y), Vector2(current_x, plot.end.y), Color("ffd166"), 1.0)
	draw_string(ThemeDB.fallback_font, plot.position + Vector2(2, 13), "%.2f" % limits.y, HORIZONTAL_ALIGNMENT_LEFT, -1, 10, Color("9eb2c1"))
	draw_string(ThemeDB.fallback_font, plot.position + Vector2(2, plot.size.y), "%.2f" % limits.x, HORIZONTAL_ALIGNMENT_LEFT, -1, 10, Color("9eb2c1"))
	if hover_index >= 0:
		_draw_hover(plot, limits, span)

func _plot_rect() -> Rect2:
	return Rect2(Vector2(9, 36), Vector2(maxf(size.x - 18.0, 1.0), maxf(size.y - 47.0, 1.0)))

func _point_for(index: int, plot: Rect2, limits: Vector2, span: float) -> Vector2:
	var safe_index := clampi(index, 0, database.rows.size() - 1)
	var x := plot.position.x + plot.size.x * float(safe_index) / float(database.rows.size() - 1)
	var y := plot.end.y - plot.size.y * float(database.rows[safe_index][column] - limits.x) / span
	return Vector2(x, y)

func _draw_hover(plot: Rect2, limits: Vector2, span: float) -> void:
	var point := _point_for(hover_index, plot, limits, span)
	var reading: PackedFloat64Array = database.rows[hover_index]
	draw_line(Vector2(point.x, plot.position.y), Vector2(point.x, plot.end.y), Color("ffffffaa"), 1.0)
	draw_circle(point, 3.5, Color.WHITE)
	var text := "t = %.3f s    %s = %.3f" % [reading[TelemetrySchema.Column.TIME], title, reading[column]]
	var text_size := ThemeDB.fallback_font.get_string_size(text, HORIZONTAL_ALIGNMENT_LEFT, -1, 12)
	var tooltip := Rect2(point + Vector2(8, -30), text_size + Vector2(12, 10))
	if tooltip.end.x > size.x - 5:
		tooltip.position.x = point.x - tooltip.size.x - 8
	if tooltip.position.y < 28:
		tooltip.position.y = point.y + 8
	draw_style_box(_tooltip_style(), tooltip)
	draw_string(ThemeDB.fallback_font, tooltip.position + Vector2(6, 15), text, HORIZONTAL_ALIGNMENT_LEFT, -1, 12, Color.WHITE)

func _panel_style() -> StyleBoxFlat:
	var style := StyleBoxFlat.new()
	style.bg_color = Color("101b25d9")
	style.border_color = Color("355166")
	style.set_border_width_all(1)
	style.corner_radius_top_left = 6
	style.corner_radius_top_right = 6
	style.corner_radius_bottom_left = 6
	style.corner_radius_bottom_right = 6
	return style

func _tooltip_style() -> StyleBoxFlat:
	var style := _panel_style()
	style.bg_color = Color("071019f5")
	style.border_color = Color("58d6ff")
	return style
