## Visual overlays attached to the supplied rocket mesh.
## Includes a visual-only RCS thruster demonstration.
extends Node3D

const TelemetrySchema = preload("res://Scripts/telemetry_schema.gd")

var cg_marker: MeshInstance3D
var cp_marker: MeshInstance3D
var thrust_arrow: MeshInstance3D
var velocity_arrow: MeshInstance3D
var fuel: MeshInstance3D
var oxidizer: MeshInstance3D
var vector_enabled := {"thrust": true, "velocity": true}

# Visual-only RCS system.
var rcs_thrusters: Array[Dictionary] = []

func _ready() -> void:
	cg_marker = _sphere(Color("56e0ff"), 0.07)
	cp_marker = _sphere(Color("ffcf5c"), 0.06)
	thrust_arrow = _cylinder(Color("ff7b38"), 0.035)
	velocity_arrow = _cylinder(Color("57e389"), 0.025)
	fuel = _tank(Color("36a8ffff"), -0.30)
	oxidizer = _tank(Color("70e090ff"), 0.28)

	_build_rcs_visuals()


func update_telemetry(row: PackedFloat64Array, initial_fuel: float, initial_ox: float) -> void:
	if row.is_empty():
		return

	# The source mesh's long axis is local +Y, so CG/CP are placed along that axis.
	cg_marker.position = Vector3(
		row[TelemetrySchema.Column.CG_Y],
		row[TelemetrySchema.Column.CG_X],
		-row[TelemetrySchema.Column.CG_Z]
	)

	cp_marker.position = Vector3(0.0, 0.75, 0.0)

	if vector_enabled["thrust"]:
		_set_arrow(
			thrust_arrow,
			Vector3(0, -0.95, 0),
			Vector3.UP,
			clampf(row[TelemetrySchema.Column.THRUST] / 15000.0, 0.10, 2.2)
		)
	else:
		thrust_arrow.visible = false

	var world_velocity := Vector3(
		row[TelemetrySchema.Column.VE],
		-row[TelemetrySchema.Column.VD],
		-row[TelemetrySchema.Column.VN]
	)

	var rocket_node := get_parent() as Node3D
	var velocity := rocket_node.global_transform.basis.inverse() * world_velocity

	if vector_enabled["velocity"]:
		_set_arrow(
			velocity_arrow,
			Vector3.ZERO,
			velocity,
			clampf(world_velocity.length() / 120.0, 0.08, 1.8)
		)
	else:
		velocity_arrow.visible = false

	_set_fill(
		fuel,
		row[TelemetrySchema.Column.FUEL_MASS],
		initial_fuel,
		-0.30
	)

	_set_fill(
		oxidizer,
		row[TelemetrySchema.Column.OX_MASS],
		initial_ox,
		0.28
	)

	# Visual-only RCS animation.
	_update_rcs(row)


func _build_rcs_visuals() -> void:
	# Four pairs of small RCS thrusters.
	#
	# Local rocket coordinates:
	#   +Y = nose
	#   -Y = engine
	#   X/Z = radial directions
	#
	# Each entry contains:
	#   position = nozzle location
	#   direction = direction the plume travels
	#
	# These are deliberately positioned outside the body so they are
	# immediately visible while we tune the visual effect.

	_add_rcs_thruster(
		Vector3(0.16, 0.35, 0.0),
		Vector3.RIGHT
	)

	_add_rcs_thruster(
		Vector3(-0.16, 0.35, 0.0),
		Vector3.LEFT
	)

	_add_rcs_thruster(
		Vector3(0.0, 0.35, 0.16),
		Vector3(0.0, 0.0, 1.0)
	)

	_add_rcs_thruster(
		Vector3(0.0, 0.35, -0.16),
		Vector3(0.0, 0.0, -1.0)
	)

	# Lower RCS cluster.
	_add_rcs_thruster(
		Vector3(0.16, -0.30, 0.0),
		Vector3.RIGHT
	)

	_add_rcs_thruster(
		Vector3(-0.16, -0.30, 0.0),
		Vector3.LEFT
	)

	_add_rcs_thruster(
		Vector3(0.0, -0.30, 0.16),
		Vector3(0.0, 0.0, 1.0)
	)

	_add_rcs_thruster(
		Vector3(0.0, -0.30, -0.16),
		Vector3(0.0, 0.0, -1.0)
	)


func _add_rcs_thruster(position: Vector3, direction: Vector3) -> void:
	var root := Node3D.new()
	root.position = position
	add_child(root)

	# Small metallic nozzle.
	var nozzle := MeshInstance3D.new()
	var nozzle_mesh := CylinderMesh.new()
	nozzle_mesh.top_radius = 0.018
	nozzle_mesh.bottom_radius = 0.045
	nozzle_mesh.height = 0.09

	nozzle.mesh = nozzle_mesh
	nozzle.material_override = _emission_material(
		Color("b8c7d1"),
		1.0
	)

	root.add_child(nozzle)

	# CylinderMesh points along local Y.
	# Rotate it so +Y points along the RCS firing direction.
	nozzle.quaternion = Quaternion(
		Vector3.UP,
		direction.normalized()
	)

	# Exhaust plume.
	var plume := MeshInstance3D.new()
	var plume_mesh := CylinderMesh.new()

	plume_mesh.top_radius = 0.0
	plume_mesh.bottom_radius = 0.025
	plume_mesh.height = 0.22

	plume.mesh = plume_mesh
	plume.material_override = _emission_material(
		Color("66d9ff"),
		0.82
	)

	plume.quaternion = Quaternion(
		Vector3.UP,
		direction.normalized()
	)

	# Move the plume away from the nozzle.
	plume.position = direction.normalized() * 0.13

	root.add_child(plume)

	# Store references so the plume can be animated.
	rcs_thrusters.append({
		"root": root,
		"nozzle": nozzle,
		"plume": plume,
		"direction": direction.normalized()
	})


func _update_rcs(row: PackedFloat64Array) -> void:
	if rcs_thrusters.is_empty():
		return

	# For now this is deliberately NOT real RCS physics.
	# We use the existing control signal to create a convincing
	# visual demonstration.
	var control := row[TelemetrySchema.Column.CONTROL]

	# Convert the control signal into an intensity.
	var intensity := clampf(absf(control), 0.0, 1.0)

	# Add a little deterministic variation so the exhaust does not
	# look completely static during playback.
	var time := row[TelemetrySchema.Column.TIME]
	var pulse := 0.72 + 0.28 * sin(time * 18.0)

	for i in rcs_thrusters.size():
		var thruster: Dictionary = rcs_thrusters[i]
		var plume: MeshInstance3D = thruster["plume"]

		# Alternate which thrusters fire based on the control direction.
		var fires := false

		if control > 0.05:
			fires = i % 2 == 0
		elif control < -0.05:
			fires = i % 2 == 1

		var strength := intensity * pulse if fires else 0.0

		plume.visible = strength > 0.04

		if plume.visible:
			var direction: Vector3 = thruster["direction"]

			plume.position = direction * (
				0.13 + strength * 0.05
			)

			plume.scale = Vector3(
				0.55 + strength * 0.45,
				0.35 + strength * 1.8,
				0.55 + strength * 0.45
			)

			var material := plume.material_override as StandardMaterial3D
			if material != null:
				material.emission_energy_multiplier = 1.5 + strength * 4.0
				material.albedo_color.a = 0.35 + strength * 0.5
		else:
			plume.scale = Vector3.ONE


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
	mesh.height = 1.0
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


func _set_arrow(
	node: MeshInstance3D,
	start: Vector3,
	direction: Vector3,
	length: float
) -> void:
	if direction.length_squared() < 0.00001:
		node.visible = false
		return

	node.visible = true
	node.position = start + direction.normalized() * length * 0.5
	node.quaternion = Quaternion(Vector3.UP, direction.normalized())
	node.scale = Vector3.ONE
	node.scale.y = length


func _set_fill(
	node: MeshInstance3D,
	mass: float,
	full_mass: float,
	local_y: float
) -> void:
	var fraction := clampf(
		mass / maxf(full_mass, 0.001),
		0.02,
		1.0
	)

	node.scale.y = fraction
	node.position.y = local_y - 0.375 * (1.0 - fraction)


func set_vector_enabled(kind: String, enabled: bool) -> void:
	if vector_enabled.has(kind):
		vector_enabled[kind] = enabled


func set_overlay_palette(palette_name: String) -> void:
	var colors := _palette_colors(palette_name)

	_set_overlay_color(cg_marker, colors["cg"])
	_set_overlay_color(cp_marker, colors["cp"])
	_set_overlay_color(thrust_arrow, colors["thrust"])
	_set_overlay_color(velocity_arrow, colors["velocity"])
	_set_overlay_color(fuel, colors["fuel"])
	_set_overlay_color(oxidizer, colors["oxidizer"])

	# Keep the RCS exhaust consistent with the overlay palette.
	var rcs_color: Color = colors["rcs"]

	for thruster in rcs_thrusters:
		var plume: MeshInstance3D = thruster["plume"]
		_set_overlay_color(plume, rcs_color)


func _palette_colors(palette_name: String) -> Dictionary:
	match palette_name:
		"high_contrast":
			return {
				"cg": Color("00e5ff"),
				"cp": Color("ffe600"),
				"thrust": Color("ff6b00"),
				"velocity": Color("65ff7a"),
				"fuel": Color("2d9cff"),
				"oxidizer": Color("65ff7a"),
				"rcs": Color("66eaff")
			}

		"monochrome":
			return {
				"cg": Color("e8eef2"),
				"cp": Color("d0d8de"),
				"thrust": Color("ffffff"),
				"velocity": Color("b9c8d1"),
				"fuel": Color("a8b8c2"),
				"oxidizer": Color("c7d3da"),
				"rcs": Color("e8eef2")
			}

		_:
			return {
				"cg": Color("56e0ff"),
				"cp": Color("ffcf5c"),
				"thrust": Color("ff7b38"),
				"velocity": Color("57e389"),
				"fuel": Color("36a8ff"),
				"oxidizer": Color("70e090"),
				"rcs": Color("66d9ff")
			}


func _set_overlay_color(node: MeshInstance3D, color: Color) -> void:
	var material := node.material_override as StandardMaterial3D

	if material == null:
		return

	material.albedo_color = Color(
		color.r,
		color.g,
		color.b,
		material.albedo_color.a
	)

	material.emission = color


func _emission_material(
	color: Color,
	alpha: float = 1.0
) -> StandardMaterial3D:
	var material := StandardMaterial3D.new()

	material.albedo_color = Color(
		color.r,
		color.g,
		color.b,
		alpha
	)

	material.emission_enabled = true
	material.emission = color
	material.emission_energy_multiplier = 1.5

	material.transparency = (
		BaseMaterial3D.TRANSPARENCY_ALPHA
		if alpha < 1.0
		else BaseMaterial3D.TRANSPARENCY_DISABLED
	)

	return material
