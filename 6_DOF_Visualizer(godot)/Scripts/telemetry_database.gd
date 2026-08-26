## In-memory telemetry history shared by playback, graphs, and the 3D view.
extends Node

const TelemetrySchema = preload("res://Scripts/telemetry_schema.gd")

#signals that can be tied into to tell other components when data is updated
signal data_changed
signal row_appended(row: PackedFloat64Array)
signal source_changed(label: String)

#number of readings before archival
const MAX_LIVE_READINGS := 1000
var rows: Array[PackedFloat64Array] = []
var source_label := "Demo trajectory"
var _tcp_server := TCPServer.new()
var _tcp_clients: Array[StreamPeerTCP] = []
var _live_enabled := false

func clear() -> void:
	rows.clear()
	data_changed.emit()

func load_csv(path: String) -> bool:
	var parsed_rows := read_csv_rows(path)
	if parsed_rows.is_empty():
		return false
	clear()
	rows = parsed_rows
	source_label = path.get_file()
	source_changed.emit(source_label)
	data_changed.emit()
	return not rows.is_empty()

func read_csv_rows(path: String) -> Array[PackedFloat64Array]:
	var file := FileAccess.open(path, FileAccess.READ)
	if file == null:
		push_error("Unable to open telemetry CSV: " + path)
		return []
	var parsed_rows: Array[PackedFloat64Array] = []
	while not file.eof_reached():
		var text := file.get_line().strip_edges()
		if not text.is_empty():
			var row := _parse_csv_tokens(text.split(",", false))
			if not row.is_empty():
				parsed_rows.append(row)
	file.close()
	return parsed_rows

func _parse_csv_tokens(tokens: PackedStringArray) -> PackedFloat64Array:
	if tokens.size() < TelemetrySchema.LEGACY_EULER_COLUMN_COUNT:
		return PackedFloat64Array()
	# Header rows cannot have a numeric timestamp.
	if not tokens[TelemetrySchema.Column.TIME].is_valid_float():
		return PackedFloat64Array()
	if tokens.size() == TelemetrySchema.LEGACY_EULER_COLUMN_COUNT:
		return _convert_legacy_euler_row(tokens)
	var result := PackedFloat64Array()
	result.resize(TelemetrySchema.COLUMN_COUNT)
	for i in TelemetrySchema.COLUMN_COUNT:
		result[i] = tokens[i].strip_edges().to_float()
	return result

func _convert_legacy_euler_row(tokens: PackedStringArray) -> PackedFloat64Array:
	var legacy_rotation := Quaternion.from_euler(Vector3(tokens[7].to_float(), tokens[5].to_float(), tokens[6].to_float()))
	var result := PackedFloat64Array()
	result.resize(TelemetrySchema.COLUMN_COUNT)
	for index in range(0, 5):
		result[index] = tokens[index].to_float()
	result[TelemetrySchema.Column.Q0] = legacy_rotation.w
	result[TelemetrySchema.Column.Q1] = legacy_rotation.x
	result[TelemetrySchema.Column.Q2] = legacy_rotation.y
	result[TelemetrySchema.Column.Q3] = legacy_rotation.z
	for index in range(9, TelemetrySchema.COLUMN_COUNT):
		result[index] = tokens[index - 1].to_float()
	return result

func append_live_csv(line: String) -> void:
	var row := _parse_csv_tokens(line.strip_edges().split(",", false))
	if not row.is_empty():
		append_live_row(row)

func append_live_row(row: PackedFloat64Array) -> void:
	if row.size() != TelemetrySchema.COLUMN_COUNT:
		push_warning("Ignored live reading with %d columns; expected %d." % [row.size(), TelemetrySchema.COLUMN_COUNT])
		return
	rows.append(row)
	if rows.size() > MAX_LIVE_READINGS:
		_archive_and_trim(rows.size() - MAX_LIVE_READINGS)
	row_appended.emit(row)
	data_changed.emit()

func start_tcp(port: int = 8080) -> bool:
	if _tcp_server.is_listening():
		return true
	var error := _tcp_server.listen(port, "*")
	_live_enabled = error == OK
	if _live_enabled:
		source_label = "Live TCP :%d" % port
		source_changed.emit(source_label)
	else:
		push_error("TCP listener could not start: %s" % error_string(error))
	return _live_enabled

func stop_tcp() -> void:
	_live_enabled = false
	_tcp_clients.clear()
	if _tcp_server.is_listening():
		_tcp_server.stop()

func _process(_delta: float) -> void:
	if not _live_enabled:
		return
	while _tcp_server.is_connection_available():
		_tcp_clients.append(_tcp_server.take_connection())
	for client in _tcp_clients.duplicate():
		if client.get_status() != StreamPeerTCP.STATUS_CONNECTED:
			_tcp_clients.erase(client)
			continue
		var count: int = client.get_available_bytes()
		if count > 0:
			var packet: String = client.get_utf8_string(count)
			for line in packet.split("\n", false):
				append_live_csv(line)

func _archive_and_trim(count: int) -> void:
	var archive_path := "user://telemetry_archive_%s.csv" % Time.get_datetime_string_from_system().replace(":", "-")
	var archive := FileAccess.open(archive_path, FileAccess.WRITE)
	if archive != null:
		archive.store_line(",".join(TelemetrySchema.FIELD_NAMES))
		for i in count:
			archive.store_line(_row_as_csv(rows[i]))
		archive.close()
	rows = rows.slice(count)

func export_csv(path: String) -> bool:
	var file := FileAccess.open(path, FileAccess.WRITE)
	if file == null:
		push_error("Unable to export telemetry to: " + path)
		return false
	file.store_line(",".join(TelemetrySchema.FIELD_NAMES))
	for row in rows:
		file.store_line(_row_as_csv(row))
	file.close()
	return true

func _row_as_csv(row: PackedFloat64Array) -> String:
	var values := PackedStringArray()
	for value in row:
		values.append(str(value))
	return ",".join(values)

func sample_at(index: int) -> PackedFloat64Array:
	return rows[clampi(index, 0, rows.size() - 1)] if not rows.is_empty() else PackedFloat64Array()

func value_range(column: int) -> Vector2:
	if rows.is_empty():
		return Vector2.ZERO
	var low := INF
	var high := -INF
	for row in rows:
		low = minf(low, row[column])
		high = maxf(high, row[column])
	return Vector2(low, high)
