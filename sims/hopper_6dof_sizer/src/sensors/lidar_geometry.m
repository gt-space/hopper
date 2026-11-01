function [d_b, r_b] = lidar_geometry(lidar)
%LIDAR_GEOMETRY  Beam directions and mounts of the four lidars, sim body axes.
%   [d_b, r_b] = lidar_geometry(S.lidar)
%     lidar  .azimuth (deg about body x, from +y), .cant (deg from the long
%            axis), .radius (m off-axis), .mount_x (m above the base)
%   d_b  3 x 4 unit beam directions (down the thrust axis, canted out)
%   r_b  3 x 4 mounts, from the base reference point
%
%   The sensor model ranges along these and the nav adapter hands the same
%   ones to the filter (nav_measurements), so the two always agree.
%
%#codegen

d_b = zeros(3, 4); r_b = zeros(3, 4);
for k = 1:4
    az = lidar.azimuth(k);
    d_b(:, k) = [-cosd(lidar.cant); sind(lidar.cant) * cosd(az); sind(lidar.cant) * sind(az)];
    r_b(:, k) = [lidar.mount_x; lidar.radius * cosd(az); lidar.radius * sind(az)];
end
end
