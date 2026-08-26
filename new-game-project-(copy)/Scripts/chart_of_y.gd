extends "res://Scripts/chart.gd"


# Called when the node enters the scene tree for the first time.
func _ready() -> void:
	var altitude_function = Function.new([], [], "Altitude", {
		type = Function.Type.LINE
	})
	self.plot([altitude_function], chart_properties)
	 # Replace with function body.


# Called every frame. 'delta' is the elapsed time since the previous frame.
func _process(delta: float) -> void:
	pass
	
func _update_graph(row):
	var time := row[0].to_float()
	var altitude := row[3].to_float()
	altitude_function.add_point(time, altitude)
	chart.queue_redraw()
