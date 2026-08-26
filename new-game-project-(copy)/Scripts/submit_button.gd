extends Button
signal csvupdated(filePath)
@onready var csvInput : TextEdit = get_node("../filePathInput")
# Called when the node enters the scene tree for the first time.
func _ready() -> void:
	print(self)
	self.pressed.connect(_on_button_pressed) # Replace with function body.


# Called every frame. 'delta' is the elapsed time since the previous frame.
func _process(delta: float) -> void:
	pass
	
func _on_button_pressed():
	var filePath: String = csvInput.text;
	print("pressed")
	csvupdated.emit(filePath)
