## Builds and updates every 3D visual element in the flight view.
extends Node

const TelemetrySchema = preload("res://Scripts/telemetry_schema.gd")
const RocketVisualScript = preload("res://Scripts/rocket_visual.gd")
const AttitudeTransform = preload("res://Scripts/attitude_transform.gd")
const ROCKET_SCENE = preload("res://Prefabs/Rocket.tscn")
const ROCKET_GROUND_OFFSET := 2.14

var rocket: Node3D
var rocket_visual: Node3D
var camera: Camera3D
var trace_mesh := ImmediateMesh.new()
var reference_mesh := ImmediateMesh.new()
var projection_mesh := ImmediateMesh.new()
var camera_mode := 0
var camera_azimuth := 35.0
var camera_elevation := 18.0

func _ready() -> void:
	_build_environment()
	_build_ground_and_rocket()
	_build_trajectory_guides()
	camera = Camera3D.new()
	camera.fov = 62.0
	add_child(camera)

func _build_environment() -> void:
	var environment := WorldEnvironment.new()
	var settings := Environment.new()
	settings.background_mode = Environment.BG_SKY
	var sky := Sky.new()
	var panorama := PanoramaSkyMaterial.new()
	panorama.panorama = load("res://grasslands_sunset_4k.hdr")
	sky.sky_material = panorama
	settings.sky = sky
	settings.ambient_light_source = Environment.AMBIENT_SOURCE_SKY
	settings.ambient_light_energy = 0.65
	settings.tonemap_mode = Environment.TONE_MAPPER_FILMIC
	environment.environment = settings
	add_child(environment)
	var sun := DirectionalLight3D.new()
	sun.rotation_degrees = Vector3(-52, -28, 0)
	sun.light_energy = 1.8
	sun.shadow_enabled = true
	add_child(sun)

func _build_ground_and_rocket() -> void:
	var ground := MeshInstance3D.new()
	var plane := PlaneMesh.new()
	plane.size = Vector2(300, 300)
	ground.mesh = plane
	ground.material_override = _material(Color("2c4932"))
	add_child(ground)
	rocket = ROCKET_SCENE.instantiate()
	rocket.scale = Vector3.ONE * 2.0
	rocket.position.y = ROCKET_GROUND_OFFSET
	add_child(rocket)
	rocket_visual = RocketVisualScript.new()
	rocket.add_child(rocket_visual)

func _build_trajectory_guides() -> void:
	_add_line_mesh(trace_mesh, Color("58d6ff"))
	_add_line_mesh(reference_mesh, Color("ffcf5c"))
	_add_line_mesh(projection_mesh, Color("5fc7ff80"))
	var ellipse := ImmediateMesh.new()
	_add_line_mesh(ellipse, Color("f2b84b"))
	ellipse.surface_begin(Mesh.PRIMITIVE_LINE_STRIP)
	for i in 65:
		var angle := TAU * float(i) / 64.0
		ellipse.surface_add_vertex(Vector3(cos(angle) * 14.0, 0.03, sin(angle) * 8.0))
	ellipse.surface_end()
	var landing_target := MeshInstance3D.new()
	var target_mesh := CylinderMesh.new()
	target_mesh.top_radius = 0.45
	target_mesh.bottom_radius = 0.45
	target_mesh.height = 0.08
	landing_target.mesh = target_mesh
	landing_target.position = Vector3(0, 0.05, 0)
	landing_target.material_override = _material(Color("ff4f64"), true)
	add_child(landing_target)
	_build_uncertainty_visuals()

func _build_uncertainty_visuals() -> void:
	var cloud := MultiMeshInstance3D.new()
	var multi_mesh := MultiMesh.new()
	multi_mesh.transform_format = MultiMesh.TRANSFORM_3D
	multi_mesh.instance_count = 60
	var dot := SphereMesh.new()
	dot.radius = 0.06
	dot.height = 0.12
	multi_mesh.mesh = dot
	for i in multi_mesh.instance_count:
		var angle := float(i) * 2.4
		var radius := 0.12 * sqrt(float(i))
		multi_mesh.set_instance_transform(i, Transform3D(Basis(), Vector3(cos(angle) * radius, 0.08, sin(angle) * radius)))
	cloud.multimesh = multi_mesh
	cloud.material_override = _material(Color("bd8cff"), true)
	cloud.name = "TrajectoryCloud"
	add_child(cloud)
	var covariance := MeshInstance3D.new()
	var ellipsoid := SphereMesh.new()
	ellipsoid.radius = 1.0
	ellipsoid.height = 2.0
	covariance.mesh = ellipsoid
	covariance.scale = Vector3(2.4, 0.35, 1.4)
	covariance.position = Vector3(0, 0.35, 0)
	covariance.material_override = _material(Color("a88cff"), true, 0.18)
	add_child(covariance)

func set_cloud_visible(visible: bool) -> void:
	get_node("TrajectoryCloud").visible = visible

func set_trace(rows: Array[PackedFloat64Array]) -> void:
	trace_mesh.clear_surfaces()
	projection_mesh.clear_surfaces()
	if rows.size() < 2:
		return
	trace_mesh.surface_begin(Mesh.PRIMITIVE_LINE_STRIP)
	projection_mesh.surface_begin(Mesh.PRIMITIVE_LINE_STRIP)
	for index in rows.size():
		var row := rows[index]
		trace_mesh.surface_add_vertex(telemetry_to_world(row))
		projection_mesh.surface_add_vertex(Vector3(row[TelemetrySchema.Column.EAST], 0.035, -row[TelemetrySchema.Column.NORTH]))
	trace_mesh.surface_end()
	projection_mesh.surface_end()

func set_reference_trace(rows: Array[PackedFloat64Array]) -> void:
	reference_mesh.clear_surfaces()
	if rows.size() < 2:
		return
	reference_mesh.surface_begin(Mesh.PRIMITIVE_LINE_STRIP)
	for row in rows:
		reference_mesh.surface_add_vertex(telemetry_to_world(row))
	reference_mesh.surface_end()

func apply_sample(_index: int, row: PackedFloat64Array, initial_fuel: float, initial_ox: float) -> void:
	rocket.position = telemetry_to_world(row)
	rocket.basis = AttitudeTransform.matlab_body_to_ned(
		row[TelemetrySchema.Column.Q0], row[TelemetrySchema.Column.Q1], row[TelemetrySchema.Column.Q2], row[TelemetrySchema.Column.Q3]
	).scaled(Vector3.ONE * 2.0)
	rocket_visual.update_telemetry(row, initial_fuel, initial_ox)

func set_vector_enabled(kind: String, enabled: bool) -> void:
	rocket_visual.set_vector_enabled(kind, enabled)

func orbit(relative_motion: Vector2) -> void:
	camera_mode = 1
	camera_azimuth -= relative_motion.x * 0.25
	camera_elevation = clampf(camera_elevation - relative_motion.y * 0.25, 5.0, 85.0)

func set_camera_mode(mode: int) -> void:
	camera_mode = mode

func update_camera(delta: float) -> void:
	var target := rocket.global_position
	var desired := target + Vector3(8, 4, 10)
	var look_target := target
	var camera_up := Vector3.UP
	if camera_mode == 1:
		var azimuth := deg_to_rad(camera_azimuth)
		var elevation := deg_to_rad(camera_elevation)
		desired = target + Vector3(cos(azimuth) * cos(elevation), sin(elevation), sin(azimuth) * cos(elevation)) * 15.0
	elif camera_mode == 2:
		desired = rocket.to_global(Vector3(0, -1.55, 0))
		desired.y = maxf(desired.y, 0.20)
		look_target = Vector3(target.x, 0.03, target.z)
		camera_up = rocket.global_transform.basis.z.normalized()
	camera.global_position = camera.global_position.lerp(desired, 1.0 - exp(-4.5 * delta))
	var direction := (look_target - camera.global_position).normalized()
	if absf(camera_up.dot(direction)) > 0.95:
		camera_up = Vector3.FORWARD
	camera.look_at(look_target, camera_up)

func telemetry_to_world(row: PackedFloat64Array) -> Vector3:
	return Vector3(row[TelemetrySchema.Column.EAST], maxf(row[TelemetrySchema.Column.UP], 0.0) + ROCKET_GROUND_OFFSET, -row[TelemetrySchema.Column.NORTH])

func _add_line_mesh(mesh: ImmediateMesh, color: Color) -> void:
	var line := MeshInstance3D.new()
	line.mesh = mesh
	var material := _material(color, true)
	material.emission_enabled = true
	material.emission = color
	material.emission_energy_multiplier = 1.2
	line.material_override = material
	add_child(line)

func _material(color: Color, unshaded: bool = false, alpha: float = 1.0) -> StandardMaterial3D:
	var material := StandardMaterial3D.new()
	material.albedo_color = Color(color, alpha)
	material.shading_mode = BaseMaterial3D.SHADING_MODE_UNSHADED if unshaded else BaseMaterial3D.SHADING_MODE_PER_PIXEL
	material.transparency = BaseMaterial3D.TRANSPARENCY_ALPHA if alpha < 1.0 else BaseMaterial3D.TRANSPARENCY_DISABLED
	return material
