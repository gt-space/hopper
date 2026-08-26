## Deterministic starter data so the visualizer can be inspected without a CSV.
extends RefCounted

static func create_rows() -> Array[PackedFloat64Array]:
	var rows: Array[PackedFloat64Array] = []
	for i in 480:
		var time := float(i) * 0.05
		# Body +X points up at launch in NED, matching the real MATLAB test flight.
		var launch_attitude := Quaternion(Vector3.UP, PI * 0.5)
		var disturbance := Quaternion.from_euler(Vector3(0.14 * sin(time * 0.45), 0.10 * sin(time), 0.07 * cos(time * 0.7)))
		var attitude := launch_attitude * disturbance
		rows.append(PackedFloat64Array([
			time,
			0.45 * time,
			0.8 * sin(time * 0.36),
			maxf(0.0, 18.0 * sin(time * 0.13)),
			10500.0 + 3600.0 * sin(time * 0.16),
			attitude.w,
			attitude.x,
			attitude.y,
			attitude.z,
			260.0 - time * 6.5,
			430.0 - time * 10.5,
			sin(time * 0.8),
			0.0,
			0.15 * sin(time),
			0.05,
			9.0,
			0.8 * cos(time * 0.36),
			-2.3 * cos(time * 0.13)
		]))
	return rows
