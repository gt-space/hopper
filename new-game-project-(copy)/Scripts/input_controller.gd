extends Node
#current line of the data stream
var currLine: int = 0
var inputReady: bool = false
var savedFilePath: String = ""
var IOType = "CSV"

#the main database for all sensor or csv data
var db = {}


#the simulation parameters
var INPUT_CONTROLLER = self
@export var simulation_paused = true
@export var timeBetweenSamples = .1
@export var stepSize: int = 5;
signal step(dbLine)

@onready var submitButton : Button = $"../UI/Canvas/CSV_Input/submitButton"
# Called when the node enters the scene tree for the first time.
func _ready() -> void:
	submitButton.csvupdated.connect(_on_CSV_updated) # Replace with function body.

#reads the csv file and converts to an array. Immediatly fills the db
func _on_CSV_updated(filePath:String):
	#pauses the sim and declares that db contains current data
	inputReady = false;
	
	#gets the filepath from the submit button
	savedFilePath = filePath;
	
	#clears and resets to beginning of sim
	db.clear()
	currLine = 0;
	
	#csv file to array
	var file = FileAccess.open(filePath, FileAccess.READ)
	if file == null:
		push_error("Could not open file")
		return null
	
	var data = []
	
	while file.get_position() < file.get_length():
		var line = file.get_line().strip_edges()
		if line == "":
			continue
		var row = line.split(",")
		#should be 18
		data.append(row)
	file.close()
	
	db = data
	#change once play button is introduced
	inputReady = true
	simulation_paused = false

# Called every frame. 'delta' is the elapsed time since the previous frame.
func _process(delta: float) -> void:
	if IOType == "CSV":
		if len(db) < currLine:
				simulation_paused = true
		if inputReady and !simulation_paused:
			step.emit(db[currLine]);
			currLine = currLine + stepSize
	#incomplete
	elif IOType == "Datastream":
		pass
#steps to a line in the database, will be used for jumping to timestamps
func _stepTo(lineNumber: int):
	if len(db) < currLine:
		simulation_paused = true
	if inputReady and !simulation_paused:
		if lineNumber > len(db) || lineNumber < 0:
			push_error("Tried to step to an index that is out of bounds")
		step.emit(db[lineNumber]);
		currLine = currLine + stepSize;
#Steps forward or backward a set number of frames from where the playback head is at TODO
func _stepForward(stepSize: int, forward: boolean):
	pass

func _togglePaused():
	simulation_paused = !simulation_paused
	
func getDB(cols: Array):
	if IOType == "CSV":
		var data = []
		for index in cols:
			data.append(db[index].slice(0, currLine))
		return data
	if IOType == "Datastream":
		pass

	
