extends Node
var currTime: int = 0
var inputReady: bool = false
var savedFilePath: String = ""
var csvFile = []
var timeBetweenSamples = .1
@export var simulation_paused = true
@export var INPUT_CONTROLLER = self
signal step(csvLine)

@onready var submitButton : Button = $"../UI/UI/submitButton"
# Called when the node enters the scene tree for the first time.
func _ready() -> void:
	submitButton.csvupdated.connect(_on_CSV_updated) # Replace with function body.
	

# Called every frame. 'delta' is the elapsed time since the previous frame.
func _process(delta: float) -> void:
	if len(csvFile) < currTime:
			simulation_paused = true
	if inputReady and !simulation_paused:
		step.emit(csvFile[currTime]);
		currTime = currTime + 5
		
func _on_CSV_updated(filePath:String):
	inputReady = true;
	savedFilePath = filePath;
	csvFile.clear()
	simulation_paused = false
	currTime = 0;
	#csv file to array
	var file = FileAccess.open(filePath, FileAccess.READ)
	print(file)
	if file == null:
		push_error("Could not open file")
		return null
	
	var data = []

	while not file.eof_reached():
		var line = file.get_line().strip_edges()
		if line == "":
			continue
		print(data)
		var row = line.split(",")
	
		data.append(row)
	csvFile = data
	
