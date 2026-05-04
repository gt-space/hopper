extends Node3D


# Called when the node enters the scene tree for the first time.
func _ready() -> void:
	get_node("../../InputController").step.connect(_on_step) # Replace with function body.


# Called every frame. 'delta' is the elapsed time since the previous frame.
func _process(delta: float) -> void:
	pass
	
func _on_step(csvLine):
	self.position.x = float(csvLine[1])
	self.position.y = float(csvLine[3])
	self.position.z = float(csvLine[2])
	
	self.rotation = Quaternion(csvLine[5], csvLine[6], csvLine[7], csvLine[8]).get_euler()
