extends Camera3D
@onready var rocket = $"../Rocket"

# Called when the node enters the scene tree for the first time.
func _ready() -> void:
	pass # Replace with function body.


# Called every frame. 'delta' is the elapsed time since the previous frame.
func _process(delta: float) -> void:
	var offset = 5
	self.position.x = rocket.position.x + offset
	self.position.y = rocket.position.y
	self.position.z = rocket.position.z
	self.look_at(rocket.global_position)
	
