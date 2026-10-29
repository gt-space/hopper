## Owns timestamp-driven playback and exposes high-level playback commands.
extends Node

signal frame_changed(start_index: int, end_index: int, interpolation: float)

const LOOP_DELAY_SECONDS := 5.0

var sample_count := 0
var index := 0
# Playback stays paused until the visualization workspace is opened.  The
# dashboard is constructed behind the Main Menu, so starting paused prevents a
# hidden mission from advancing before the user presses Load Mission.
var is_playing := false
var speed := 1.0
var loop_enabled := false
var _timestamps := PackedFloat64Array()
var _playhead := 0.0
var _loop_delay_remaining := -1.0

func set_rows(rows: Array[PackedFloat64Array]) -> void:
	var previous_count := sample_count
	sample_count = rows.size()
	_timestamps.resize(sample_count)
	for row_index in sample_count:
		var timestamp := rows[row_index][0]
		# A zero-duration or out-of-order sample cannot be interpolated safely.
		if row_index > 0:
			timestamp = maxf(timestamp, _timestamps[row_index - 1] + 0.000001)
		_timestamps[row_index] = timestamp
	index = clampi(index, 0, maxi(sample_count - 1, 0))
	if previous_count == 0 and sample_count > 0:
		_playhead = _timestamps[0]
	elif sample_count > 0:
		_playhead = clampf(_playhead, _timestamps[0], _timestamps[-1])

func reset() -> void:
	restart_playing()


func rewind_and_pause() -> void:
	_rewind()
	is_playing = false
	emit_current_frame()


func restart_playing() -> void:
	_rewind()
	is_playing = true
	emit_current_frame()


func _rewind() -> void:
	index = 0
	_playhead = _timestamps[0] if not _timestamps.is_empty() else 0.0
	_loop_delay_remaining = -1.0

func toggle_playing() -> void:
	is_playing = not is_playing

func set_speed(multiplier: float) -> void:
	speed = multiplier

func set_loop_enabled(enabled: bool) -> void:
	loop_enabled = enabled
	if not loop_enabled:
		_loop_delay_remaining = -1.0

func step(amount: int) -> void:
	is_playing = false
	_loop_delay_remaining = -1.0
	index = clampi(index + amount, 0, maxi(sample_count - 1, 0))
	if not _timestamps.is_empty():
		_playhead = _timestamps[index]
	emit_current_frame()

func scrub_to(value: float) -> void:
	_loop_delay_remaining = -1.0
	index = clampi(roundi(value), 0, maxi(sample_count - 1, 0))
	if not _timestamps.is_empty():
		_playhead = _timestamps[index]
	emit_current_frame()

func emit_current_frame() -> void:
	if sample_count > 0:
		frame_changed.emit(index, index, 0.0)

func _process(delta: float) -> void:
	if not is_playing or sample_count < 2:
		return
	if _loop_delay_remaining >= 0.0:
		_loop_delay_remaining -= delta
		if _loop_delay_remaining <= 0.0:
			index = 0
			_playhead = _timestamps[0]
			_loop_delay_remaining = -1.0
			emit_current_frame()
		return
	_playhead += delta * speed
	if _playhead >= _timestamps[-1]:
		_playhead = _timestamps[-1]
		index = sample_count - 1
		emit_current_frame()
		if loop_enabled:
			_loop_delay_remaining = LOOP_DELAY_SECONDS
		else:
			is_playing = false
		return
	# Advancing directly to the current timestamp intentionally skips intermediate
	# display frames at high playback speeds instead of doing work for each sample.
	while index < sample_count - 2 and _timestamps[index + 1] <= _playhead:
		index += 1
	var next_index := index + 1
	var interval := maxf(_timestamps[next_index] - _timestamps[index], 0.000001)
	var interpolation := clampf((_playhead - _timestamps[index]) / interval, 0.0, 1.0)
	frame_changed.emit(index, next_index, interpolation)
