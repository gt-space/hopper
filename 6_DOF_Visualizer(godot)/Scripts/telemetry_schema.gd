## Canonical row layout for MATLAB/Simulink telemetry.
extends RefCounted

## MATLAB's A, B, C, D quaternion export is assumed scalar-first: [q0, q1, q2, q3].
enum Column { TIME, NORTH, EAST, UP, THRUST, Q0, Q1, Q2, Q3, FUEL_MASS, OX_MASS, CONTROL, CG_X, CG_Y, CG_Z, VN, VE, VD }

const COLUMN_COUNT := 18
const LEGACY_EULER_COLUMN_COUNT := 17
const FIELD_NAMES = ["time", "north", "east", "up", "thrust", "q0", "q1", "q2", "q3", "fuel_mass", "ox_mass", "control", "cg_x", "cg_y", "cg_z", "vn", "ve", "vd"]
const DISPLAY_NAMES = ["Time", "North", "East", "Altitude Up", "Thrust", "Quaternion q0", "Quaternion q1", "Quaternion q2", "Quaternion q3", "Fuel mass", "Oxidizer mass", "Control", "CG X", "CG Y", "CG Z", "Velocity North", "Velocity East", "Velocity Down"]

static func display_name(index: int) -> String:
	return DISPLAY_NAMES[index] if index >= 0 and index < DISPLAY_NAMES.size() else "Unknown"
