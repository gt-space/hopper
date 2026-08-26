## Converts scalar-first MATLAB body-to-NED quaternions into Godot's native basis.
extends RefCounted

## NED -> native Godot: (north, east, down) -> (east, up, south) = (east, -down, -north).
## This is a proper right-handed mapping and leaves Godot's normal XZ ground plane intact.
const NED_TO_GODOT := Basis(Vector3(0, 0, -1), Vector3(1, 0, 0), Vector3(0, -1, 0))
## Supplied mesh local axes: +Y is the rocket nose. MATLAB body axes: +X forward, +Y right, +Z down.
const BODY_FROM_MODEL := Basis(Vector3(0, 1, 0), Vector3(1, 0, 0), Vector3(0, 0, -1))

static func matlab_body_to_ned(q0: float, q1: float, q2: float, q3: float) -> Basis:
	var quaternion := Quaternion(q1, q2, q3, q0).normalized()
	return NED_TO_GODOT * Basis(quaternion) * BODY_FROM_MODEL
