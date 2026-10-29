## Owns persistent, runtime-loadable vehicle profiles.  Custom files never alter res://.
extends Node

signal selected_profile_changed(profile: Dictionary)
signal profiles_changed

const VEHICLES_PATH := "user://vehicles"
const SELECTION_PATH := "user://vehicles/selection.cfg"
const DEFAULT_PROFILE := {
	"id": "hopper-default",
	"display_name": "Hopper — Default Vehicle",
	"created_unix": 0,
	"model_path": "",
	"dimensions": {"height": 4.88, "width": 1.23, "length": 1.23},
	"rcs_offset": {"x": 0.0, "y": 0.0, "z": 0.0},
	"tvc_offset": {"x": 0.0, "y": 0.0, "z": 0.0},
	"model_scale": 2.0,
	"model_rotation": {"x": 0.0, "y": 0.0, "z": 0.0},
}

var profiles: Array[Dictionary] = []
var selected_id := "hopper-default"

func load_from_disk() -> void:
	DirAccess.make_dir_recursive_absolute(ProjectSettings.globalize_path(VEHICLES_PATH))
	profiles = [DEFAULT_PROFILE.duplicate(true)]
	var directory := DirAccess.open(VEHICLES_PATH)
	if directory != null:
		for folder in directory.get_directories():
			var profile := _read_profile("%s/%s/profile.json" % [VEHICLES_PATH, folder])
			if not profile.is_empty() and FileAccess.file_exists(profile.get("model_path", "")):
				profiles.append(profile)
	profiles.sort_custom(func(a: Dictionary, b: Dictionary): return int(a["created_unix"]) < int(b["created_unix"]))
	var selection := ConfigFile.new()
	if selection.load(SELECTION_PATH) == OK:
		selected_id = str(selection.get_value("vehicle", "selected_id", "hopper-default"))
	if get_selected_profile().is_empty():
		selected_id = "hopper-default"
		_save_selection()
	profiles_changed.emit()
	selected_profile_changed.emit(get_selected_profile())

func get_selected_profile() -> Dictionary:
	return get_profile(selected_id)

func get_profile(id: String) -> Dictionary:
	for profile in profiles:
		if profile.get("id", "") == id:
			return profile.duplicate(true)
	return {}

func select_profile(id: String) -> bool:
	var profile := get_profile(id)
	if profile.is_empty():
		return false
	selected_id = id
	_save_selection()
	selected_profile_changed.emit(profile)
	return true

func save_profile(profile: Dictionary) -> bool:
	var id := str(profile.get("id", ""))
	if id.is_empty() or id == DEFAULT_PROFILE["id"]:
		return false
	var path := "%s/%s/profile.json" % [VEHICLES_PATH, id]
	var file := FileAccess.open(path, FileAccess.WRITE)
	if file == null:
		return false
	file.store_string(JSON.stringify(profile, "\t"))
	file.close()
	for index in profiles.size():
		if profiles[index].get("id", "") == id:
			profiles[index] = profile.duplicate(true)
			break
	profiles_changed.emit()
	if id == selected_id:
		selected_profile_changed.emit(profile.duplicate(true))
	return true

func create_profile_from_glb(source_path: String) -> Dictionary:
	if source_path.get_extension().to_lower() != "glb" or not FileAccess.file_exists(source_path):
		return {}
	var probe := _load_glb(source_path)
	if probe == null:
		return {}
	var dimensions := _estimate_dimensions(probe)
	probe.queue_free()
	var stem := source_path.get_file().get_basename().validate_filename()
	if stem.is_empty():
		stem = "vehicle"
	var id := "%s_%d" % [stem.to_lower(), int(Time.get_unix_time_from_system())]
	var profile_dir := "%s/%s" % [VEHICLES_PATH, id]
	var suffix := 2
	while DirAccess.dir_exists_absolute(ProjectSettings.globalize_path(profile_dir)):
		id = "%s_%d" % [id, suffix]
		profile_dir = "%s/%s" % [VEHICLES_PATH, id]
		suffix += 1
	if DirAccess.make_dir_recursive_absolute(ProjectSettings.globalize_path(profile_dir)) != OK:
		return {}
	var target_path := "%s/model.glb" % profile_dir
	if DirAccess.copy_absolute(source_path, ProjectSettings.globalize_path(target_path)) != OK:
		return {}
	var profile := {
		"id": id,
		"display_name": stem,
		"created_unix": int(Time.get_unix_time_from_system()),
		"model_path": target_path,
		"dimensions": dimensions,
		"rcs_offset": {"x": 0.0, "y": 0.0, "z": 0.0},
		"tvc_offset": {"x": 0.0, "y": 0.0, "z": 0.0},
		"model_scale": 1.0,
		"model_rotation": {"x": 0.0, "y": 0.0, "z": 0.0},
	}
	if not save_profile(profile):
		return {}
	profiles.append(profile.duplicate(true))
	profiles.sort_custom(func(a: Dictionary, b: Dictionary): return int(a["created_unix"]) < int(b["created_unix"]))
	profiles_changed.emit()
	select_profile(id)
	return profile

func instantiate_model(profile: Dictionary) -> Node3D:
	var model: Node3D
	if str(profile.get("id", "")) == DEFAULT_PROFILE["id"]:
		model = preload("res://Prefabs/Rocket.tscn").instantiate()
	else:
		model = _load_glb(str(profile.get("model_path", "")))
	if model == null:
		return null
	# Imported CAD exports often have geometry offset from their root.  Put the
	# mesh beneath a stable pivot and center its lateral bounds on that pivot so
	# carousel rotation and flight overlays share the intended vehicle axis.
	var bounds := _mesh_bounds(model)
	var center: Vector3 = (bounds["minimum"] as Vector3 + bounds["maximum"] as Vector3) * 0.5
	model.position -= Vector3(center.x, 0.0, center.z)
	var pivot := Node3D.new()
	pivot.name = "VehicleModelPivot"
	pivot.add_child(model)
	pivot.scale = Vector3.ONE * float(profile.get("model_scale", 1.0))
	var rotation: Dictionary = profile.get("model_rotation", {})
	pivot.rotation_degrees = Vector3(float(rotation.get("x", 0.0)), float(rotation.get("y", 0.0)), float(rotation.get("z", 0.0)))
	return pivot

func _load_glb(path: String) -> Node3D:
	var document := GLTFDocument.new()
	var state := GLTFState.new()
	if document.append_from_file(path, state) != OK:
		return null
	var generated := document.generate_scene(state)
	if generated is Node3D:
		return generated as Node3D
	if generated != null:
		generated.queue_free()
	return null

func _estimate_dimensions(root: Node) -> Dictionary:
	var bounds := _mesh_bounds(root)
	var size: Vector3 = bounds["maximum"] - bounds["minimum"]
	return {"height": size.y, "width": size.x, "length": size.z}

func _mesh_bounds(root: Node) -> Dictionary:
	var bounds := {
		"minimum": Vector3(1000000.0, 1000000.0, 1000000.0),
		"maximum": Vector3(-1000000.0, -1000000.0, -1000000.0),
		"found": false,
	}
	_collect_mesh_bounds(root, Transform3D.IDENTITY, bounds)
	if not bounds["found"]:
		bounds["minimum"] = Vector3.ZERO
		bounds["maximum"] = Vector3.ZERO
	return bounds

func _collect_mesh_bounds(node: Node, parent_transform: Transform3D, bounds: Dictionary) -> void:
	var node_transform := parent_transform
	if node is Node3D:
		node_transform = parent_transform * (node as Node3D).transform
	if node is MeshInstance3D:
		var mesh := (node as MeshInstance3D).mesh
		if mesh != null:
			var aabb := mesh.get_aabb()
			for x in [aabb.position.x, aabb.end.x]:
				for y in [aabb.position.y, aabb.end.y]:
					for z in [aabb.position.z, aabb.end.z]:
						var point := node_transform * Vector3(x, y, z)
						var minimum: Vector3 = bounds["minimum"]
						var maximum: Vector3 = bounds["maximum"]
						bounds["minimum"] = Vector3(minf(minimum.x, point.x), minf(minimum.y, point.y), minf(minimum.z, point.z))
						bounds["maximum"] = Vector3(maxf(maximum.x, point.x), maxf(maximum.y, point.y), maxf(maximum.z, point.z))
						bounds["found"] = true
	for child in node.get_children():
		_collect_mesh_bounds(child, node_transform, bounds)

func _read_profile(path: String) -> Dictionary:
	var file := FileAccess.open(path, FileAccess.READ)
	if file == null:
		return {}
	var json := JSON.new()
	if json.parse(file.get_as_text()) != OK or not json.data is Dictionary:
		return {}
	return json.data as Dictionary

func _save_selection() -> void:
	var config := ConfigFile.new()
	config.set_value("vehicle", "selected_id", selected_id)
	config.save(SELECTION_PATH)
