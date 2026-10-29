function names = fdir_targets()
% FDIR_TARGETS  Task targets a fault's response acts on (fault table column).
%
% R0 masks go out as a bitmask in this order (bit k-1 for target k), so
% navigation and the allocator can drop the source.
names = ["imu", "imu_accel", "imu_gyro", "lidar1", "lidar2", "lidar3", "lidar4", "lidar", ...
         "gps", "baro", "mag", "nav", "guidance", "control", "tvc", "rcs", "throttle"];
end
