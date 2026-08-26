## Visual overlays attached to the supplied rocket mesh.
extends Node3D

const TelemetrySchema = preload("res://Scripts/telemetry_schema.gd")

var cg_marker: MeshInstance3D
var cp_marker: MeshInstance3D
var thrust_arrow: MeshInstance3D
var velocity_arrow: MeshInstance3D
var rcs_arrow: MeshInstance3D
var fuel: MeshInstance3D
var oxidizer: MeshInstance3D
var vector_enabled := {"thrust": true, "velocity": true, "rcs": true}

func _ready() -> void:
	cg_marker = _sphere(Color("56e0ff"), 0.07)
	cp_marker = _sphere(Color("ffcf5c"), 0.06)
	thrust_arrow = _cylinder(Color("ff7b38"), 0.035)
	velocity_arrow = _cylinder(Color("57e389"), 0.025)
	rcs_arrow = _cylinder(Color("e85aad"), 0.02)
	fuel = _tank(Color("36a8ffff"), -0.30)
	oxidizer = _tank(Color("70e090ff"), 0.28)

func update_telemetry(row: PackedFloat64Array, initial_fuel: float, initial_ox: float) -> void:
	if row.is_empty():
		return
	# The source mesh's long axis is local +Y, so CG/CP are placed along that axis.
	cg_marker.position = Vector3(row[TelemetrySchema.Column.CG_Y], row[TelemetrySchema.Column.CG_X], -row[TelemetrySchema.Column.CG_Z])
	cp_marker.position = Vector3(0.0, 0.75, 0.0)
	if vector_enabled["thrust"]:
		_set_arrow(thrust_arrow, Vector3(0, -0.95, 0), Vector3.UP, clampf(row[TelemetrySchema.Column.THRUST] / 15000.0, 0.10, 2.2))
	else:
		thrust_arrow.visible = false
	var world_velocity := Vector3(row[TelemetrySchema.Column.VE], -row[TelemetrySchema.Column.VD], -row[TelemetrySchema.Column.VN])
	var rocket_node := get_parent() as Node3D
	var velocity := rocket_node.global_transform.basis.inverse() * world_velocity
	if vector_enabled["velocity"]:
		_set_arrow(velocity_arrow, Vector3.ZERO, velocity, clampf(world_velocity.length() / 120.0, 0.08, 1.8))
	else:
		velocity_arrow.visible = false
	var command := row[TelemetrySchema.Column.CONTROL]
	if vector_enabled["rcs"]:
		_set_arrow(rcs_arrow, Vector3(0.22, 0.3, 0), Vector3(command, 0.15, 0), clampf(absf(command), 0.04, 0.8))
	else:
		rcs_arrow.visible = false
	_set_fill(fuel, row[TelemetrySchema.Column.FUEL_MASS], initial_fuel, -0.30)
	_set_fill(oxidizer, row[TelemetrySchema.Column.OX_MASS], initial_ox, 0.28)

func _sphere(color: Color, radius: float) -> MeshInstance3D:
	var node := MeshInstance3D.new()
	var mesh := SphereMesh.new()
	mesh.radius = radius
	mesh.height = radius * 2.0
	node.mesh = mesh
	node.material_override = _emission_material(color)
	add_child(node)
	return node

func _cylinder(color: Color, radius: float) -> MeshInstance3D:
	var node := MeshInstance3D.new()
	var mesh := CylinderMesh.new()
	mesh.top_radius = radius * 0.3
	mesh.bottom_radius = radius
	node.mesh = mesh
	node.material_override = _emission_material(color)
	add_child(node)
	return node

func _tank(color: Color, local_y: float) -> MeshInstance3D:
	var node := MeshInstance3D.new()
	var mesh := CylinderMesh.new()
	mesh.top_radius = 0.095
	mesh.bottom_radius = 0.095
	mesh.height = 0.75
	node.mesh = mesh
	node.material_override = _emission_material(color, 0.45)
	node.position = Vector3(0, local_y, 0)
	add_child(node)
	return node

func _set_arrow(node: MeshInstance3D, start: Vector3, direction: Vector3, length: float) -> void:
	if direction.length_squared() < 0.00001:
		node.visible = false
		return
	node.visible = true
	node.position = start + direction.normalized() * length * 0.5
	node.quaternion = Quaternion(Vector3.UP, direction.normalized())
	node.scale = Vector3.ONE
	node.scale.y = length

func _set_fill(node: MeshInstance3D, mass: float, full_mass: float, local_y: float) -> void:
	var fraction := clampf(mass / maxf(full_mass, 0.001), 0.02, 1.0)
	node.scale.y = fraction
	node.position.y = local_y - 0.375 * (1.0 - fraction)

func set_vector_enabled(kind: String, enabled: bool) -> void:
	if vector_enabled.has(kind):
		vector_enabled[kind] = enabled

func _emission_material(color: Color, alpha: float = 1.0) -> StandardMaterial3D:
	var material := StandardMaterial3D.new()
	material.albedo_color = Color(color, alpha)
	material.emission_enabled = true
	material.emission = color
	material.emission_energy_multiplier = 1.5
	material.transparency = BaseMaterial3D.TRANSPARENCY_ALPHA if alpha < 1.0 else BaseMaterial3D.TRANSPARENCY_DISABLED
	return material
