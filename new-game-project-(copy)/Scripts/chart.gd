extends Chart

@onready var InputController = get_node("root/InputController")
var chart = self
var cp
var f1
# Called when the node enters the scene tree for the first time.
func _ready() -> void:
	InputController.step.connect(_refresh())
	var x: Array = []
	
	# And our y values. It can be an n-size array of arrays.
	# NOTE: `x.size() == y.size()` or `x.size() == y[n].size()`
	var y: Array = []
	
	# Let's customize the chart properties, which specify how the chart
	# should look, plus some additional elements like labels, the scale, etc...
	cp = ChartProperties.new()
	cp.colors.frame = Color("#161a1d")
	cp.colors.background = Color.TRANSPARENT
	cp.colors.grid = Color("#283442")
	cp.colors.ticks = Color("#283442")
	cp.colors.text = Color.WHITE_SMOKE
	cp.y_scale = 10
	cp.draw_origin = true
	cp.draw_bounding_box = false
	cp.draw_vertical_grid = false
	cp.interactive = true # false by default, it allows the chart to create a tooltip to show point values
	# and interecept clicks on the plot
	
	# Let's add values to our functions
	f1 = Function.new(
		x, y, "User", # This will create a function with x and y values taken by the Arrays 
						# we have created previously. This function will also be named "Pressure"
						# as it contains 'pressure' values.
						# If set, the name of a function will be used both in the Legend
						# (if enabled thourgh ChartProperties) and on the Tooltip (if enabled).
		{
			type = Function.Type.BAR,
			bar_size = 5
		}
	)
	
	# Now let's plot our data
	chart.plot([f1], cp)
	
	# Uncommenting this line will show how real time data plotting works
	set_process(false)


# Called every frame. 'delta' is the elapsed time since the previous frame.
func _process(delta: float) -> void:
	pass

func _refresh():
	var data = InputController.getDB([0,3])
	chart.datasets = [{ "x": data[0], "y": data[1] }]
	chart.plot([f1], cp)
	
