## Builds and updates every 3D visual element in the flight view.
extends Node

const TelemetrySchema = preload("res://Scripts/telemetry_schema.gd")
const RocketVisualScript = preload("res://Scripts/rocket_visual.gd")
const EnvironmentTerrainScript = preload("res://Scripts/environment_terrain.gd")
const AttitudeTransform = preload("res://Scripts/attitude_transform.gd")
const ROCKET_SCENE = preload("res://Prefabs/Rocket.tscn")
const ROCKET_GROUND_OFFSET := 2.14
const LANDING_CAMERA_LOCAL_OFFSET := Vector3(0.0, -0.92, 0.0)
const LANDING_CAMERA_LOCAL_BASIS := Basis(Vector3.RIGHT, Vector3.FORWARD, Vector3.UP)

var rocket: Node3D
var rocket_visual: Node3D
var vehicle_model: Node3D
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
var environment_terrain: Node
var landing_ellipse_mesh := ImmediateMesh.new()
var camera_azimuth := 35.0
var camera_elevation := 18.0
var _camera_states: Dictionary = {}
var _ground_origin_y := 0.0
var _trace_rows: Array[PackedFloat64Array] = []
var _reference_rows: Array[PackedFloat64Array] = []
var _current_row := PackedFloat64Array()

func _ready() -> void:
	_build_environment()
	_build_ground_and_rocket()
	_build_trajectory_guides()
	environment_terrain = EnvironmentTerrainScript.new()
	environment_terrain.status_changed.connect(func(message: String): environment_status_changed.emit(message))
	environment_terrain.ground_origin_changed.connect(_on_ground_origin_changed)
	add_child(environment_terrain)

signal environment_status_changed(message: String)

func _build_environment() -> void:
	var environment := WorldEnvironment.new()
	var settings := Environment.new()
	settings.background_mode = Environment.BG_SKY
	var sky := Sky.new()
	var sky_material := ProceduralSkyMaterial.new()
	sky_material.sky_top_color = Color("287fb7")
	sky_material.sky_horizon_color = Color("8fc8ce")
	sky_material.ground_bottom_color = Color("1a322a")
	sky_material.ground_horizon_color = Color("5f9878")
	sky_material.sun_angle_max = 9.0
	sky_material.sun_curve = 0.05
	sky.sky_material = sky_material
	settings.sky = sky
	settings.ambient_light_source = Environment.AMBIENT_SOURCE_SKY
	settings.ambient_light_energy = 0.54
	settings.fog_enabled = true
	settings.fog_light_color = Color("9fc0ae")
	settings.fog_density = 0.00085
	settings.fog_sky_affect = 0.32
	settings.tonemap_mode = Environment.TONE_MAPPER_FILMIC
	environment.environment = settings
	add_child(environment)
	var sun := DirectionalLight3D.new()
	sun.rotation_degrees = Vector3(-52, -28, 0)
	sun.light_energy = 2.15
	sun.shadow_enabled = true
	add_child(sun)

func _build_ground_and_rocket() -> void:
	# Keep the flight transform on a stable pivot so a profile model can be
	# replaced without touching its telemetry position, attitude, or overlays.
	rocket = Node3D.new()
	rocket.name = "VehiclePivot"
	rocket.position.y = ROCKET_GROUND_OFFSET
	add_child(rocket)
	_set_vehicle_model(ROCKET_SCENE.instantiate())

func set_vehicle_model(model: Node3D) -> void:
	if model == null:
		return
	_set_vehicle_model(model)

func _set_vehicle_model(model: Node3D) -> void:
	if is_instance_valid(vehicle_model):
		rocket.remove_child(vehicle_model)
		vehicle_model.queue_free()
	vehicle_model = model
	rocket.add_child(vehicle_model)
	rocket_visual = RocketVisualScript.new()
	vehicle_model.add_child(rocket_visual)

func _build_trajectory_guides() -> void:
	trace_line = _add_line_mesh(trace_mesh, Color("58d6ff"))
	reference_line = _add_line_mesh(reference_mesh, Color("ffcf5c"))
	projection_line = _add_line_mesh(projection_mesh, Color("ff334f"))
	landing_ellipse = _add_line_mesh(landing_ellipse_mesh, Color("f2b84b"))
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
	_refresh_ground_references()

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
	set_visual_colors(_palette_colors(palette_name))

func set_visual_colors(colors: Dictionary) -> void:
	_set_overlay_color(trace_line, colors["trace"])
	_set_overlay_color(reference_line, colors["reference"])
	_set_overlay_color(projection_line, colors["projection"])
	_set_overlay_color(landing_ellipse, colors["landing"])
	_set_overlay_color(landing_target, colors["landing"])
	_set_overlay_color(ground_position_marker, colors["projection"])
	_set_overlay_color(monte_carlo_cloud, colors["uncertainty"])
	_set_overlay_color(covariance_ellipsoid, colors["uncertainty"])
	if is_instance_valid(rocket_visual):
		rocket_visual.set_overlay_colors(colors)

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
	_trace_rows = rows.duplicate()
	trace_mesh.clear_surfaces()
	projection_mesh.clear_surfaces()
	if rows.size() < 2:
		return
	trace_mesh.surface_begin(Mesh.PRIMITIVE_TRIANGLES)
	projection_mesh.surface_begin(Mesh.PRIMITIVE_LINE_STRIP)
	for index in range(1, rows.size()):
		_append_thin_ribbon_segment(trace_mesh, telemetry_to_world(rows[index - 1]), telemetry_to_world(rows[index]))
	for row in rows:
		var x := row[TelemetrySchema.Column.EAST]
		var z := -row[TelemetrySchema.Column.NORTH]
		projection_mesh.surface_add_vertex(Vector3(x, _surface_height_at(x, z) + 0.035, z))
	trace_mesh.surface_end()
	projection_mesh.surface_end()

func set_reference_trace(rows: Array[PackedFloat64Array]) -> void:
	_reference_rows = rows.duplicate()
	reference_mesh.clear_surfaces()
	if rows.size() < 2:
		return
	reference_mesh.surface_begin(Mesh.PRIMITIVE_LINE_STRIP)
	for row in rows:
		reference_mesh.surface_add_vertex(telemetry_to_world(row))
	reference_mesh.surface_end()

func apply_sample(_index: int, row: PackedFloat64Array, initial_fuel: float, initial_ox: float) -> void:
	_current_row = row
	rocket.position = telemetry_to_world(row)
	rocket.basis = AttitudeTransform.matlab_body_to_ned(
		row[TelemetrySchema.Column.Q0], row[TelemetrySchema.Column.Q1], row[TelemetrySchema.Column.Q2], row[TelemetrySchema.Column.Q3]
	)
	ground_position_marker.position = Vector3(rocket.position.x, _surface_height_at(rocket.position.x, rocket.position.z) + 0.14, rocket.position.z)
	rocket_visual.update_telemetry(row, initial_fuel, initial_ox)

func set_vector_enabled(kind: String, enabled: bool) -> void:
	rocket_visual.set_vector_enabled(kind, enabled)

func load_terrain(latitude: float, longitude: float, api_key: String, area_size_km: float) -> void:
	if is_instance_valid(environment_terrain):
		environment_terrain.request_live_terrain(latitude, longitude, api_key, area_size_km)

func restore_default_terrain() -> void:
	if is_instance_valid(environment_terrain):
		environment_terrain.restore_default()

func orbit(relative_motion: Vector2) -> void:
	camera_azimuth -= relative_motion.x * 0.25
	camera_elevation = clampf(camera_elevation - relative_motion.y * 0.25, 5.0, 85.0)

func update_view_cameras(delta: float, cameras: Dictionary) -> void:
	for mode in cameras:
		var view_camera := cameras[mode] as Camera3D
		if is_instance_valid(view_camera):
			_update_view_camera(int(mode), view_camera, delta)

func _update_view_camera(mode: int, view_camera: Camera3D, delta: float) -> void:
	if mode == 2:
		_snap_landing_camera(view_camera)
		return
	var target := rocket.global_position
	var desired := target + Vector3(8, 4, 10)
	var look_target := target
	var camera_up := Vector3.UP
	if mode == 1:
		var azimuth := deg_to_rad(camera_azimuth)
		var elevation := deg_to_rad(camera_elevation)
		desired = target + Vector3(cos(azimuth) * cos(elevation), sin(elevation), sin(azimuth) * cos(elevation)) * 15.0
	var smoothing := 1.0 - exp(-4.5 * delta)
	var state: Dictionary = _camera_states.get(mode, {})
	if state.is_empty():
		view_camera.global_position = desired
		state = {"look_target": look_target, "camera_up": camera_up}
	else:
		view_camera.global_position = view_camera.global_position.lerp(desired, smoothing)
		state["look_target"] = (state["look_target"] as Vector3).lerp(look_target, smoothing)
		state["camera_up"] = (state["camera_up"] as Vector3).lerp(camera_up, smoothing).normalized()
	var smoothed_target := state["look_target"] as Vector3
	var smoothed_up := state["camera_up"] as Vector3
	var direction := (smoothed_target - view_camera.global_position).normalized()
	if absf(smoothed_up.dot(direction)) > 0.95:
		smoothed_up = Vector3.FORWARD
		state["camera_up"] = smoothed_up
	view_camera.look_at(smoothed_target, smoothed_up)
	_camera_states[mode] = state

func _snap_landing_camera(view_camera: Camera3D) -> void:
	# This view is rigidly mounted near the engine: no world-space clamp or smoothing.
	var mounted_transform := rocket.global_transform * Transform3D(LANDING_CAMERA_LOCAL_BASIS, LANDING_CAMERA_LOCAL_OFFSET)
	view_camera.global_transform = mounted_transform.orthonormalized()

func _on_ground_origin_changed(height: float) -> void:
	_ground_origin_y = height
	_refresh_ground_references()
	if not _trace_rows.is_empty():
		set_trace(_trace_rows)
	if not _reference_rows.is_empty():
		set_reference_trace(_reference_rows)
	if _current_row.is_empty():
		rocket.position.y = _ground_origin_y + ROCKET_GROUND_OFFSET
	else:
		rocket.position = telemetry_to_world(_current_row)
		ground_position_marker.position = Vector3(rocket.position.x, _surface_height_at(rocket.position.x, rocket.position.z) + 0.14, rocket.position.z)

func _refresh_ground_references() -> void:
	landing_ellipse_mesh.clear_surfaces()
	landing_ellipse_mesh.surface_begin(Mesh.PRIMITIVE_LINE_STRIP)
	for i in 65:
		var angle := TAU * float(i) / 64.0
		var x := cos(angle) * 14.0
		var z := sin(angle) * 8.0
		landing_ellipse_mesh.surface_add_vertex(Vector3(x, _surface_height_at(x, z) + 0.03, z))
	landing_ellipse_mesh.surface_end()
	landing_target.position = Vector3(0.0, _surface_height_at(0.0, 0.0) + 0.05, 0.0)
	ground_position_marker.position = Vector3(rocket.position.x, _surface_height_at(rocket.position.x, rocket.position.z) + 0.14, rocket.position.z)
	covariance_ellipsoid.position = Vector3(0.0, _surface_height_at(0.0, 0.0) + 0.35, 0.0)

func _surface_height_at(x: float, z: float) -> float:
	if is_instance_valid(environment_terrain):
		return environment_terrain.get_surface_height_at(x, z)
	return _ground_origin_y

func telemetry_to_world(row: PackedFloat64Array) -> Vector3:
	return Vector3(row[TelemetrySchema.Column.EAST], maxf(row[TelemetrySchema.Column.UP], 0.0) + ROCKET_GROUND_OFFSET + _ground_origin_y, -row[TelemetrySchema.Column.NORTH])

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
