extends Node

@onready var inputController = $InputController
@onready var rocket = $Objects/Rocket
@onready var altitudeChart = $UI/Canvas/Right_Panel/AltitudeChart
# Called when the node enters the scene tree for the first time.
func _ready() -> void:
	#connects the input controller to the rocket step function
	inputController.step.connect(rocket._on_step)
	inputController.step.connnect(altitudeChart._update_graph)

# Called every frame. 'delta' is the elapsed time since the previous frame.
func _process(delta: float) -> void:
	pass
