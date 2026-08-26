extends Node3D


# Called when the node enters the scene tree for the first time.
func _ready() -> void:
	get_node("/root/InputController").step.connect(_on_step) # Replace with function body.


# Called every frame. 'delta' is the elapsed time since the previous frame.
func _process(delta: float) -> void:
	pass
	
func _on_step(csvLine):
	self.position.x = float(csvLine[1])
	self.position.y = float(csvLine[3])
	self.position.z = float(csvLine[2])

	var rotations = Quaternion(float(csvLine[5]), float(csvLine[6]), float(csvLine[7]), float(csvLine[8])).get_euler()
	print(rotations)
	#fix this at somepoint
	self.rotation = Vector3(rotations.x,-rotations.z,rotations.y - PI/2)
