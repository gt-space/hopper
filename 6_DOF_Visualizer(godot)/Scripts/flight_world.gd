## Builds and updates every 3D visual element in the flight view.
extends Node

const TelemetrySchema = preload("res://Scripts/telemetry_schema.gd")
const RocketVisualScript = preload("res://Scripts/rocket_visual.gd")
const AttitudeTransform = preload("res://Scripts/attitude_transform.gd")
const ROCKET_SCENE = preload("res://Prefabs/Rocket.tscn")
const ROCKET_GROUND_OFFSET := 2.14
const LANDING_CAMERA_LOCAL_OFFSET := Vector3(0.0, -0.92, 0.0)
const LANDING_CAMERA_LOCAL_BASIS := Basis(Vector3.RIGHT, Vector3.FORWARD, Vector3.UP)

var rocket: Node3D
var rocket_visual: Node3D
var camera: Camera3D
var trace_mesh := ImmediateMesh.new()
var reference_mesh := ImmediateMesh.new()
var projection_mesh := ImmediateMesh.new()
var trace_line: MeshInstance3D
var reference_line: MeshInstance3D
var projection_line: MeshInstance3D
var ground_position_marker: MeshInstance3D
var landing_ellipse: MeshInstance3D
var landing_target: MeshInstance3D
var monte_carlo_cloud: MultiMeshInstance3D
var covariance_ellipsoid: MeshInstance3D
var camera_mode := 0
var camera_azimuth := 35.0
var camera_elevation := 18.0
var _smoothed_look_target := Vector3.ZERO
var _smoothed_camera_up := Vector3.UP
var _camera_initialized := false

func _ready() -> void:
	_build_environment()
	_build_ground_and_rocket()
	_build_trajectory_guides()
	_build_cloud_layer()
	camera = Camera3D.new()
	camera.fov = 62.0
	add_child(camera)

func _build_environment() -> void:
	var environment := WorldEnvironment.new()
	var settings := Environment.new()
	settings.background_mode = Environment.BG_SKY
	var sky := Sky.new()
	var panorama := PanoramaSkyMaterial.new()
	panorama.panorama = load("res://Assets/Skies/cartoon_cloud_sky.png")
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

func _build_cloud_layer() -> void:
	var cloud_layer := Node3D.new()
	cloud_layer.name = "CloudLayer100m"
	cloud_layer.position.y = 100.0
	add_child(cloud_layer)
	var puff_mesh := SphereMesh.new()
	puff_mesh.radius = 1.0
	puff_mesh.height = 2.0
	var cloud_material := _material(Color("f7fbff"), true, 0.18)
	cloud_material.cull_mode = BaseMaterial3D.CULL_DISABLED
	cloud_material.emission_enabled = true
	cloud_material.emission = Color("d8eeff")
	cloud_material.emission_energy_multiplier = 0.35
	for i in 48:
		var puff := MeshInstance3D.new()
		puff.mesh = puff_mesh
		puff.material_override = cloud_material
		var angle := float(i) * 2.39996
		var radius := 15.0 + sqrt(float(i)) * 18.0
		puff.position = Vector3(cos(angle) * radius, sin(float(i) * 1.7) * 0.35, sin(angle) * radius)
		var stretch := 4.5 + fmod(float(i) * 1.73, 3.5)
		puff.scale = Vector3(stretch, 0.16 + fmod(float(i), 3.0) * 0.025, stretch * 0.65)
		cloud_layer.add_child(puff)

func _build_trajectory_guides() -> void:
	trace_line = _add_line_mesh(trace_mesh, Color("58d6ff"))
	reference_line = _add_line_mesh(reference_mesh, Color("ffcf5c"))
	projection_line = _add_line_mesh(projection_mesh, Color("ff334f"))
	var ellipse := ImmediateMesh.new()
	landing_ellipse = _add_line_mesh(ellipse, Color("f2b84b"))
	ellipse.surface_begin(Mesh.PRIMITIVE_LINE_STRIP)
	for i in 65:
		var angle := TAU * float(i) / 64.0
		ellipse.surface_add_vertex(Vector3(cos(angle) * 14.0, 0.03, sin(angle) * 8.0))
	ellipse.surface_end()
	landing_target = MeshInstance3D.new()
	var target_mesh := CylinderMesh.new()
	target_mesh.top_radius = 0.45
	target_mesh.bottom_radius = 0.45
	target_mesh.height = 0.08
	landing_target.mesh = target_mesh
	landing_target.position = Vector3(0, 0.05, 0)
	landing_target.material_override = _material(Color("ff4f64"), true)
	add_child(landing_target)
	ground_position_marker = MeshInstance3D.new()
	var marker_mesh := SphereMesh.new()
	marker_mesh.radius = 0.13
	marker_mesh.height = 0.26
	ground_position_marker.mesh = marker_mesh
	ground_position_marker.position.y = 0.14
	ground_position_marker.material_override = _material(Color("ff334f"), true)
	var marker_material := ground_position_marker.material_override as StandardMaterial3D
	marker_material.emission_enabled = true
	marker_material.emission = Color("ff334f")
	marker_material.emission_energy_multiplier = 2.5
	add_child(ground_position_marker)
	_build_uncertainty_visuals()

func _build_uncertainty_visuals() -> void:
	monte_carlo_cloud = MultiMeshInstance3D.new()
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
	monte_carlo_cloud.multimesh = multi_mesh
	monte_carlo_cloud.material_override = _material(Color("bd8cff"), true)
	add_child(monte_carlo_cloud)
	covariance_ellipsoid = MeshInstance3D.new()
	var ellipsoid := SphereMesh.new()
	ellipsoid.radius = 1.0
	ellipsoid.height = 2.0
	covariance_ellipsoid.mesh = ellipsoid
	covariance_ellipsoid.scale = Vector3(2.4, 0.35, 1.4)
	covariance_ellipsoid.position = Vector3(0, 0.35, 0)
	covariance_ellipsoid.material_override = _material(Color("a88cff"), true, 0.18)
	add_child(covariance_ellipsoid)

func set_landing_ellipse_visible(visible: bool) -> void:
	landing_ellipse.visible = visible
	landing_target.visible = visible

func set_monte_carlo_visible(visible: bool) -> void:
	monte_carlo_cloud.visible = visible
	covariance_ellipsoid.visible = visible

func set_overlay_palette(palette_name: String) -> void:
	var colors := _palette_colors(palette_name)
	_set_overlay_color(trace_line, colors["trace"])
	_set_overlay_color(reference_line, colors["reference"])
	_set_overlay_color(projection_line, colors["projection"])
	_set_overlay_color(landing_ellipse, colors["landing"])
	_set_overlay_color(landing_target, colors["landing"])
	_set_overlay_color(ground_position_marker, colors["projection"])
	_set_overlay_color(monte_carlo_cloud, colors["uncertainty"])
	_set_overlay_color(covariance_ellipsoid, colors["uncertainty"])
	rocket_visual.set_overlay_palette(palette_name)

func _palette_colors(palette_name: String) -> Dictionary:
	match palette_name:
		"high_contrast":
			return {"trace": Color("00e5ff"), "reference": Color("ffe600"), "projection": Color("ff2b63"), "landing": Color("ff8c00"), "uncertainty": Color("b15cff")}
		"monochrome":
			return {"trace": Color("d7e4ec"), "reference": Color("b6c5cf"), "projection": Color("ffffff"), "landing": Color("e6e6e6"), "uncertainty": Color("96a6b2")}
		_:
			return {"trace": Color("58d6ff"), "reference": Color("ffcf5c"), "projection": Color("ff334f"), "landing": Color("f2b84b"), "uncertainty": Color("bd8cff")}

func _set_overlay_color(node: GeometryInstance3D, color: Color) -> void:
	var material := node.material_override as StandardMaterial3D
	if material == null:
		return
	material.albedo_color = Color(color.r, color.g, color.b, material.albedo_color.a)
	if material.emission_enabled:
		material.emission = color

func set_trace(rows: Array[PackedFloat64Array]) -> void:
	trace_mesh.clear_surfaces()
	projection_mesh.clear_surfaces()
	if rows.size() < 2:
		return
	trace_mesh.surface_begin(Mesh.PRIMITIVE_TRIANGLES)
	projection_mesh.surface_begin(Mesh.PRIMITIVE_LINE_STRIP)
	for index in range(1, rows.size()):
		_append_thin_ribbon_segment(trace_mesh, telemetry_to_world(rows[index - 1]), telemetry_to_world(rows[index]))
	for row in rows:
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
	ground_position_marker.position = Vector3(rocket.position.x, 0.14, rocket.position.z)
	rocket_visual.update_telemetry(row, initial_fuel, initial_ox)
	if camera_mode == 2:
		_snap_landing_camera()

func set_vector_enabled(kind: String, enabled: bool) -> void:
	rocket_visual.set_vector_enabled(kind, enabled)

func orbit(relative_motion: Vector2) -> void:
	camera_mode = 1
	camera_azimuth -= relative_motion.x * 0.25
	camera_elevation = clampf(camera_elevation - relative_motion.y * 0.25, 5.0, 85.0)

func set_camera_mode(mode: int) -> void:
	camera_mode = mode
	if mode == 2:
		_snap_landing_camera()
	else:
		_camera_initialized = false

func update_camera(delta: float) -> void:
	if camera_mode == 2:
		_snap_landing_camera()
		return
	var target := rocket.global_position
	var desired := target + Vector3(8, 4, 10)
	var look_target := target
	var camera_up := Vector3.UP
	if camera_mode == 1:
		var azimuth := deg_to_rad(camera_azimuth)
		var elevation := deg_to_rad(camera_elevation)
		desired = target + Vector3(cos(azimuth) * cos(elevation), sin(elevation), sin(azimuth) * cos(elevation)) * 15.0
	var smoothing := 1.0 - exp(-4.5 * delta)
	if not _camera_initialized:
		camera.global_position = desired
		_smoothed_look_target = look_target
		_smoothed_camera_up = camera_up
		_camera_initialized = true
	else:
		camera.global_position = camera.global_position.lerp(desired, smoothing)
		_smoothed_look_target = _smoothed_look_target.lerp(look_target, smoothing)
		_smoothed_camera_up = _smoothed_camera_up.lerp(camera_up, smoothing).normalized()
	var direction := (_smoothed_look_target - camera.global_position).normalized()
	if absf(_smoothed_camera_up.dot(direction)) > 0.95:
		_smoothed_camera_up = Vector3.FORWARD
	camera.look_at(_smoothed_look_target, _smoothed_camera_up)

func _snap_landing_camera() -> void:
	if camera == null:
		return
	# This view is rigidly mounted near the engine: no world-space clamp or smoothing.
	var mounted_transform := rocket.global_transform * Transform3D(LANDING_CAMERA_LOCAL_BASIS, LANDING_CAMERA_LOCAL_OFFSET)
	camera.global_transform = mounted_transform.orthonormalized()
	_camera_initialized = false

func telemetry_to_world(row: PackedFloat64Array) -> Vector3:
	return Vector3(row[TelemetrySchema.Column.EAST], maxf(row[TelemetrySchema.Column.UP], 0.0) + ROCKET_GROUND_OFFSET, -row[TelemetrySchema.Column.NORTH])

func _append_thin_ribbon_segment(mesh: ImmediateMesh, start: Vector3, finish: Vector3) -> void:
	var direction := finish - start
	if direction.length_squared() < 0.000001:
		return
	var side := direction.cross(Vector3.UP)
	if side.length_squared() < 0.000001:
		side = Vector3.RIGHT
	side = side.normalized() * 0.018
	mesh.surface_add_vertex(start - side)
	mesh.surface_add_vertex(start + side)
	mesh.surface_add_vertex(finish + side)
	mesh.surface_add_vertex(start - side)
	mesh.surface_add_vertex(finish + side)
	mesh.surface_add_vertex(finish - side)

func _add_line_mesh(mesh: ImmediateMesh, color: Color) -> MeshInstance3D:
	var line := MeshInstance3D.new()
	line.mesh = mesh
	var material := _material(color, true)
	material.emission_enabled = true
	material.emission = color
	material.emission_energy_multiplier = 1.2
	material.cull_mode = BaseMaterial3D.CULL_DISABLED
	line.material_override = material
	add_child(line)
	return line

func _material(color: Color, unshaded: bool = false, alpha: float = 1.0) -> StandardMaterial3D:
	var material := StandardMaterial3D.new()
	material.albedo_color = Color(color, alpha)
	material.shading_mode = BaseMaterial3D.SHADING_MODE_UNSHADED if unshaded else BaseMaterial3D.SHADING_MODE_PER_PIXEL
	material.transparency = BaseMaterial3D.TRANSPARENCY_ALPHA if alpha < 1.0 else BaseMaterial3D.TRANSPARENCY_DISABLED
	return material
## Updates the landing ellipse geometry dynamically based on target or covariance statistics.
func set_landing_ellipse_params(center_east: float, center_north: float, semi_major: float, semi_minor: float, rotation_rad: float = 0.0, segments: int = 64) -> void:
	var mesh := landing_ellipse.mesh as ImmediateMesh
	if mesh == null:
		return
		
	mesh.clear_surfaces()
	mesh.surface_begin(Mesh.PRIMITIVE_LINE_STRIP)
	
	for i in range(segments + 1):
		var angle := TAU * float(i) / float(segments)
		
		# Parametric equation of an ellipse in local coordinates
		var local_x := cos(angle) * semi_major
		var local_z := sin(angle) * semi_minor
		
		# Rotate the ellipse if a rotation angle is provided
		var rot_x := local_x * cos(rotation_rad) - local_z * sin(rotation_rad)
		var rot_z := local_x * sin(rotation_rad) + local_z * cos(rotation_rad)
		
		# Map to Godot coordinates (East -> X, North -> -Z)
		var world_x := center_east + rot_x
		var world_z := -(center_north + rot_z)
		
		mesh.surface_add_vertex(Vector3(world_x, 0.04, world_z))
		
	mesh.surface_end()
