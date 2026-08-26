## Owns playback timing and exposes high-level playback commands.
extends Node

signal sample_changed(index: int)

const SAMPLE_INTERVAL := 0.05

var sample_count := 0
var index := 0
var is_playing := true
var speed := 1.0
var _elapsed := 0.0

func set_sample_count(count: int) -> void:
	sample_count = maxi(count, 0)
	index = clampi(index, 0, maxi(sample_count - 1, 0))
	_elapsed = 0.0

func reset() -> void:
	index = 0
	_elapsed = 0.0
	is_playing = true
	emit_current_sample()

func toggle_playing() -> void:
	is_playing = not is_playing

func set_speed(multiplier: float) -> void:
	speed = multiplier

func step(amount: int) -> void:
	is_playing = false
	index = clampi(index + amount, 0, maxi(sample_count - 1, 0))
	emit_current_sample()

func scrub_to(value: float) -> void:
	index = clampi(roundi(value), 0, maxi(sample_count - 1, 0))
	emit_current_sample()

func emit_current_sample() -> void:
	if sample_count > 0:
		sample_changed.emit(index)

func _process(delta: float) -> void:
	if not is_playing or sample_count < 2:
		return
	_elapsed += delta * speed
	while _elapsed >= SAMPLE_INTERVAL:
		index = mini(index + 1, sample_count - 1)
		_elapsed -= SAMPLE_INTERVAL
		if index == sample_count - 1:
			is_playing = false
		emit_current_sample()
