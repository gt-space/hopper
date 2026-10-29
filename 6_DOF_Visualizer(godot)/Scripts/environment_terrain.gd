## Owns the decorative terrain layer. It never changes flight telemetry coordinates.
extends Node3D

signal status_changed(message: String)
signal ground_origin_changed(height: float)

const GRASS_TEXTURE = preload("res://Assets/Environment/stylized_grassland_tile.png")
const DEFAULT_RANGE_SIZE := 300.0
const MAX_GRID_SAMPLES := 96
const MIN_LIVE_RANGE_KM := 1.0
const MAX_LIVE_RANGE_KM := 5.0

var _environment_root: Node3D
var _http_request: HTTPRequest
var _request_in_flight := false
var _pending_range_size := DEFAULT_RANGE_SIZE
var _active_columns := 0
var _active_rows := 0
var _active_heights := PackedFloat32Array()
var _active_range_size := DEFAULT_RANGE_SIZE

func _ready() -> void:
	_http_request = HTTPRequest.new()
	_http_request.timeout = 25.0
	add_child(_http_request)
	_http_request.request_completed.connect(_on_terrain_request_completed)
	build_test_range()

func build_test_range(announce_ready: bool = true) -> void:
	_reset_environment_root()
	var heights := PackedFloat32Array()
	var samples := 49
	for row in samples:
		for column in samples:
			var x := lerpf(-DEFAULT_RANGE_SIZE * 0.5, DEFAULT_RANGE_SIZE * 0.5, float(column) / float(samples - 1))
			var z := lerpf(-DEFAULT_RANGE_SIZE * 0.5, DEFAULT_RANGE_SIZE * 0.5, float(row) / float(samples - 1))
			heights.append(_test_range_height(x, z))
	_set_active_terrain(samples, samples, heights, DEFAULT_RANGE_SIZE)
	_build_terrain_mesh(samples, samples, heights, DEFAULT_RANGE_SIZE, "TestRangeTerrain")
	_build_launch_pad()
	_build_range_grid()
	_build_range_markers()
	_build_perimeter_props()
	_build_cloud_layer()
	ground_origin_changed.emit(_terrain_height_at(0.0, 0.0))
	if announce_ready:
		status_changed.emit("Stylized test range ready.")

func request_live_terrain(latitude: float, longitude: float, api_key: String, area_size_km: float) -> void:
	if _request_in_flight:
		status_changed.emit("Terrain request already in progress.")
		return
	if api_key.strip_edges().is_empty():
		status_changed.emit("Enter an OpenTopography API key to load terrain.")
		return
	if latitude < -90.0 or latitude > 90.0 or longitude < -180.0 or longitude > 180.0:
		status_changed.emit("Terrain coordinates are outside the valid latitude/longitude range.")
		return
	var range_size := clampf(area_size_km, MIN_LIVE_RANGE_KM, MAX_LIVE_RANGE_KM) * 1000.0
	var half_span_lat := (range_size * 0.5) / 111320.0
	var longitude_scale := maxf(cos(deg_to_rad(latitude)), 0.01)
	var half_span_lon := (range_size * 0.5) / (111320.0 * longitude_scale)
	var url := "https://portal.opentopography.org/API/globaldem?demtype=COP30&south=%.8f&north=%.8f&west=%.8f&east=%.8f&outputFormat=AAIGrid&API_Key=%s" % [latitude - half_span_lat, latitude + half_span_lat, longitude - half_span_lon, longitude + half_span_lon, api_key.uri_encode()]
	_pending_range_size = range_size
	_request_in_flight = true
	var request_error := _http_request.request(url)
	if request_error != OK:
		_request_in_flight = false
		build_test_range(false)
		status_changed.emit("Could not start terrain download (%s); restored the test range." % error_string(request_error))
		return
	status_changed.emit("Downloading %.0f × %.0f km terrain data…" % [range_size / 1000.0, range_size / 1000.0])

func restore_default() -> void:
	if _request_in_flight:
		_http_request.cancel_request()
		_request_in_flight = false
	build_test_range()

func _on_terrain_request_completed(result: int, response_code: int, _headers: PackedStringArray, body: PackedByteArray) -> void:
	_request_in_flight = false
	if result != HTTPRequest.RESULT_SUCCESS or response_code != 200:
		build_test_range(false)
		var reason := "HTTP %d" % response_code if response_code > 0 else "network result %d" % result
		status_changed.emit("Terrain download was unavailable (%s); restored the test range." % reason)
		return
	var parsed := _parse_arc_ascii_grid(body.get_string_from_utf8())
	if parsed.is_empty():
		build_test_range(false)
		status_changed.emit("Terrain data could not be read; restored the test range.")
		return
	_reset_environment_root()
	_set_active_terrain(parsed["columns"], parsed["rows"], parsed["heights"], _pending_range_size)
	_build_terrain_mesh(parsed["columns"], parsed["rows"], parsed["heights"], _pending_range_size, "LiveTerrain")
	_build_launch_pad()
	_build_range_grid()
	_build_range_markers()
	_build_cloud_layer()
	ground_origin_changed.emit(_terrain_height_at(0.0, 0.0))
	status_changed.emit("Loaded %.0f × %.0f km OpenTopography terrain." % [_pending_range_size / 1000.0, _pending_range_size / 1000.0])

func _parse_arc_ascii_grid(text: String) -> Dictionary:
	var lines := text.split("\n", false)
	if lines.size() < 7:
		return {}
	var header := {}
	for index in 6:
		var tokens := lines[index].strip_edges().split(" ", false)
		if tokens.size() < 2:
			return {}
		header[tokens[0].to_lower()] = tokens[tokens.size() - 1].to_float()
	var source_columns := int(header.get("ncols", 0.0))
	var source_rows := int(header.get("nrows", 0.0))
	if source_columns < 2 or source_rows < 2:
		return {}
	var no_data: float = float(header.get("nodata_value", -9999.0))
	var values := PackedFloat32Array()
	for index in range(6, lines.size()):
		for token in lines[index].strip_edges().replace("\t", " ").split(" ", false):
			if token.is_valid_float():
				values.append(token.to_float())
	# Some Global DEM responses advertise one more row than their AAIGrid body
	# contains.  Use only complete received rows rather than rejecting a valid
	# terrain request because of that harmless header mismatch.
	var received_rows := values.size() / source_columns
	if received_rows < 2:
		return {}
	source_rows = mini(source_rows, received_rows)
	var stride := maxi(1, ceili(float(maxi(source_columns, source_rows)) / float(MAX_GRID_SAMPLES)))
	var columns := maxi(2, ceili(float(source_columns - 1) / float(stride)) + 1)
	var rows := maxi(2, ceili(float(source_rows - 1) / float(stride)) + 1)
	var center_row := source_rows / 2
	var center_column := source_columns / 2
	var center_value := values[center_row * source_columns + center_column]
	if is_equal_approx(center_value, no_data):
		return {}
	var heights := PackedFloat32Array()
	for row in rows:
		var source_row := mini(row * stride, source_rows - 1)
		for column in columns:
			var source_column := mini(column * stride, source_columns - 1)
			var value := values[source_row * source_columns + source_column]
			# Preserve the DEM's physical vertical scale: one Godot unit is one meter.
			# Heights are still relative to the launch-point sample, not sea level.
			heights.append(0.0 if is_equal_approx(value, no_data) else value - center_value)
	return {"columns": columns, "rows": rows, "heights": heights}

func _reset_environment_root() -> void:
	if is_instance_valid(_environment_root):
		_environment_root.free()
	_environment_root = Node3D.new()
	_environment_root.name = "EnvironmentTerrain"
	add_child(_environment_root)

func _set_active_terrain(columns: int, rows: int, heights: PackedFloat32Array, range_size: float) -> void:
	_active_columns = columns
	_active_rows = rows
	_active_heights = heights
	_active_range_size = range_size

func _build_terrain_mesh(columns: int, rows: int, heights: PackedFloat32Array, range_size: float, node_name: String) -> void:
	var vertices := PackedVector3Array()
	var normals := PackedVector3Array()
	var uvs := PackedVector2Array()
	var indices := PackedInt32Array()
	var step_x := range_size / float(columns - 1)
	var step_z := range_size / float(rows - 1)
	for row in rows:
		for column in columns:
			var index := row * columns + column
			vertices.append(Vector3(-range_size * 0.5 + float(column) * step_x, heights[index], -range_size * 0.5 + float(row) * step_z))
			uvs.append(Vector2(float(column) / float(columns - 1), float(row) / float(rows - 1)))
			var left := heights[row * columns + maxi(column - 1, 0)]
			var right := heights[row * columns + mini(column + 1, columns - 1)]
			var north := heights[maxi(row - 1, 0) * columns + column]
			var south := heights[mini(row + 1, rows - 1) * columns + column]
			normals.append(Vector3(-(right - left) / maxf(step_x * 2.0, 0.001), 1.0, -(south - north) / maxf(step_z * 2.0, 0.001)).normalized())
	for row in range(rows - 1):
		for column in range(columns - 1):
			var top_left := row * columns + column
			var top_right := top_left + 1
			var bottom_left := (row + 1) * columns + column
			var bottom_right := bottom_left + 1
			indices.append_array(PackedInt32Array([top_left, bottom_left, top_right, top_right, bottom_left, bottom_right]))
	var arrays := []
	arrays.resize(ArrayMesh.ARRAY_MAX)
	arrays[ArrayMesh.ARRAY_VERTEX] = vertices
	arrays[ArrayMesh.ARRAY_NORMAL] = normals
	arrays[ArrayMesh.ARRAY_TEX_UV] = uvs
	arrays[ArrayMesh.ARRAY_INDEX] = indices
	var mesh := ArrayMesh.new()
	mesh.add_surface_from_arrays(Mesh.PRIMITIVE_TRIANGLES, arrays)
	var terrain := MeshInstance3D.new()
	terrain.name = node_name
	terrain.mesh = mesh
	terrain.material_override = _terrain_material(range_size / 18.0)
	_environment_root.add_child(terrain)

func _terrain_material(repetitions: float) -> StandardMaterial3D:
	var material := StandardMaterial3D.new()
	material.albedo_texture = GRASS_TEXTURE
	material.roughness = 0.82
	material.uv1_scale = Vector3(repetitions, repetitions, 1.0)
	material.albedo_color = Color.WHITE
	return material

func _build_launch_pad() -> void:
	var pad := MeshInstance3D.new()
	pad.name = "LaunchPad"
	var mesh := CylinderMesh.new()
	mesh.top_radius = 3.5
	mesh.bottom_radius = 3.5
	mesh.height = 0.14
	pad.mesh = mesh
	pad.position.y = _terrain_height_at(0.0, 0.0) + 0.07
	pad.material_override = _solid_material(Color("728479"))
	_environment_root.add_child(pad)
	var ring := MeshInstance3D.new()
	var ring_mesh := TorusMesh.new()
	ring_mesh.inner_radius = 3.55
	ring_mesh.outer_radius = 3.72
	ring.mesh = ring_mesh
	ring.position.y = _terrain_height_at(0.0, 0.0) + 0.10
	ring.material_override = _solid_material(Color("f0c35b"), true)
	_environment_root.add_child(ring)

func _build_range_grid() -> void:
	var grid := ImmediateMesh.new()
	grid.surface_begin(Mesh.PRIMITIVE_LINES)
	var half_range := int(_active_range_size * 0.5)
	var major_spacing := 25 if is_equal_approx(_active_range_size, DEFAULT_RANGE_SIZE) else 50
	# A closer sample spacing keeps the rendered line locked to steep DEM facets
	# after switching to a true 1:1 vertical scale.
	var line_sample_spacing := 10 if is_equal_approx(_active_range_size, DEFAULT_RANGE_SIZE) else 5
	for coordinate in range(-half_range, half_range + 1, major_spacing):
		for segment in range(-half_range, half_range, line_sample_spacing):
			var x0 := float(segment)
			var x1 := float(segment + line_sample_spacing)
			var z0 := float(segment)
			var z1 := float(segment + line_sample_spacing)
			grid.surface_add_vertex(Vector3(float(coordinate), _terrain_height_at(float(coordinate), z0) + 0.12, z0))
			grid.surface_add_vertex(Vector3(float(coordinate), _terrain_height_at(float(coordinate), z1) + 0.12, z1))
			grid.surface_add_vertex(Vector3(x0, _terrain_height_at(x0, float(coordinate)) + 0.12, float(coordinate)))
			grid.surface_add_vertex(Vector3(x1, _terrain_height_at(x1, float(coordinate)) + 0.12, float(coordinate)))
	grid.surface_end()
	var grid_instance := MeshInstance3D.new()
	grid_instance.name = "RangeGrid"
	grid_instance.mesh = grid
	grid_instance.material_override = _solid_material(Color("dae5c5"), true, 0.34)
	_environment_root.add_child(grid_instance)

func _build_range_markers() -> void:
	var marker_material := _solid_material(Color("f5bf3f"), true)
	marker_material.emission_enabled = true
	marker_material.emission = Color("b86b12")
	marker_material.emission_energy_multiplier = 0.25
	var base_material := _solid_material(Color("16344b"), true)
	var marker_directions := [
		{"direction": Vector3.RIGHT, "cardinal": "E"},
		{"direction": Vector3.LEFT, "cardinal": "W"},
		{"direction": Vector3.FORWARD, "cardinal": "S"},
		{"direction": Vector3.BACK, "cardinal": "N"},
	]
	for distance in [50.0, 100.0]:
		for marker_info in marker_directions:
			var direction: Vector3 = marker_info.direction
			var ground_position: Vector3 = direction * distance
			var ground_height := _terrain_height_at(ground_position.x, ground_position.z)

			var base := MeshInstance3D.new()
			var base_mesh := CylinderMesh.new()
			base_mesh.top_radius = 0.72
			base_mesh.bottom_radius = 0.86
			base_mesh.height = 0.12
			base.mesh = base_mesh
			base.position = ground_position + Vector3.UP * (ground_height + 0.06)
			base.material_override = base_material
			_environment_root.add_child(base)

			var marker := MeshInstance3D.new()
			var pole := CylinderMesh.new()
			pole.top_radius = 0.15
			pole.bottom_radius = 0.19
			pole.height = 3.6
			marker.mesh = pole
			marker.position = ground_position + Vector3.UP * (ground_height + pole.height * 0.5)
			marker.material_override = marker_material
			_environment_root.add_child(marker)
			var label := Label3D.new()
			label.text = "%s  %d m" % [marker_info.cardinal, int(distance)]
			label.position = marker.position + Vector3.UP * 2.15
			label.billboard = BaseMaterial3D.BILLBOARD_ENABLED
			# Use a world-sized label rather than a screen-fixed one.  It grows
			# naturally as the camera approaches, while its roughly 2 m maximum
			# physical height stays comparable to the small perimeter trees.
			label.fixed_size = false
			label.font_size = 52
			label.pixel_size = 0.04
			label.outline_size = 7
			label.outline_modulate = Color("102337")
			label.modulate = Color("fff5ca")
			_environment_root.add_child(label)

func _build_perimeter_props() -> void:
	var rng := RandomNumberGenerator.new()
	rng.seed = 6012026
	var trunk_mesh := CylinderMesh.new()
	trunk_mesh.top_radius = 0.12
	trunk_mesh.bottom_radius = 0.18
	trunk_mesh.height = 1.6
	var canopy_mesh := CylinderMesh.new()
	canopy_mesh.top_radius = 0.0
	canopy_mesh.bottom_radius = 1.25
	canopy_mesh.height = 3.8
	for index in 34:
		var angle := rng.randf_range(0.0, TAU)
		var radius := rng.randf_range(104.0, 142.0)
		var position := Vector3(cos(angle) * radius, 0.0, sin(angle) * radius)
		var scale := rng.randf_range(0.65, 1.35)
		var trunk := MeshInstance3D.new()
		trunk.mesh = trunk_mesh
		trunk.position = position + Vector3.UP * 0.8 * scale
		trunk.scale = Vector3(scale, scale, scale)
		trunk.material_override = _solid_material(Color("74513b"))
		_environment_root.add_child(trunk)
		var canopy := MeshInstance3D.new()
		canopy.mesh = canopy_mesh
		canopy.position = position + Vector3.UP * 3.0 * scale
		canopy.scale = Vector3(scale, scale, scale)
		canopy.material_override = _solid_material(Color("2e704d"))
		_environment_root.add_child(canopy)

func _build_cloud_layer() -> void:
	var clouds := Node3D.new()
	clouds.name = "CloudLayer100m"
	clouds.position.y = 100.0
	_environment_root.add_child(clouds)
	var puff_mesh := SphereMesh.new()
	puff_mesh.radius = 1.0
	puff_mesh.height = 2.0
	var material := _solid_material(Color("eaf5ff"), true, 0.72)
	material.cull_mode = BaseMaterial3D.CULL_DISABLED
	material.emission_enabled = true
	material.emission = Color("b8dcf2")
	material.emission_energy_multiplier = 0.08
	for index in 36:
		var puff := MeshInstance3D.new()
		puff.mesh = puff_mesh
		puff.material_override = material
		var angle := float(index) * 2.39996
		var radius := 46.0 + sqrt(float(index)) * 18.0
		puff.position = Vector3(cos(angle) * radius, sin(float(index) * 1.7) * 0.65, sin(angle) * radius)
		var stretch := 4.2 + fmod(float(index) * 1.73, 3.0)
		puff.scale = Vector3(stretch, 0.30 + fmod(float(index), 3.0) * 0.04, stretch * 0.74)
		clouds.add_child(puff)

func _terrain_height_at(x: float, z: float) -> float:
	if _active_columns < 2 or _active_rows < 2 or _active_heights.is_empty():
		return 0.0
	var normalized_x := clampf((x + _active_range_size * 0.5) / _active_range_size, 0.0, 1.0) * float(_active_columns - 1)
	var normalized_z := clampf((z + _active_range_size * 0.5) / _active_range_size, 0.0, 1.0) * float(_active_rows - 1)
	var left := int(floor(normalized_x))
	var top := int(floor(normalized_z))
	var right := mini(left + 1, _active_columns - 1)
	var bottom := mini(top + 1, _active_rows - 1)
	var blend_x := normalized_x - float(left)
	var blend_z := normalized_z - float(top)
	var top_left := _active_heights[top * _active_columns + left]
	var top_right := _active_heights[top * _active_columns + right]
	var bottom_left := _active_heights[bottom * _active_columns + left]
	var bottom_right := _active_heights[bottom * _active_columns + right]
	# Match the exact two-triangle winding used in _build_terrain_mesh().  A
	# bilinear sample does not lie on a triangulated mesh except at its vertices.
	if blend_x + blend_z <= 1.0:
		return top_left + blend_x * (top_right - top_left) + blend_z * (bottom_left - top_left)
	return bottom_right + (1.0 - blend_z) * (top_right - bottom_right) + (1.0 - blend_x) * (bottom_left - bottom_right)

func get_surface_height_at(x: float, z: float) -> float:
	return _terrain_height_at(x, z)

func _test_range_height(x: float, z: float) -> float:
	# The offline range is intentionally level so the launch pad, reference grid,
	# and telemetry origin all read as a single, unambiguous ground plane.
	return 0.0

func _solid_material(color: Color, unshaded: bool = false, alpha: float = 1.0) -> StandardMaterial3D:
	var material := StandardMaterial3D.new()
	material.albedo_color = Color(color.r, color.g, color.b, alpha)
	material.roughness = 0.88
	material.shading_mode = BaseMaterial3D.SHADING_MODE_UNSHADED if unshaded else BaseMaterial3D.SHADING_MODE_PER_PIXEL
	material.transparency = BaseMaterial3D.TRANSPARENCY_ALPHA if alpha < 1.0 else BaseMaterial3D.TRANSPARENCY_DISABLED
	return material
