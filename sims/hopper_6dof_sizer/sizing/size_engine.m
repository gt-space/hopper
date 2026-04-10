function ENGINE = size_engine(IN)
%% Description
% 1) Unpack Inputs - nominal thrust, throttle range, Pe, C* efficiency
% size for fully expanded at 20% stiffness at 100% thrust

% 2) Use imported engine contour and calculate TCA mass with
% regen_thickness - multiply by 1.2 to account for manifold

% 3) Use chamber diameter from engine contour to calculate injector volume
% and mass

% 4) Use thrust, gimbal range, and gimbal velocity to calculate TVC
% actuator mass - multiply by 2X

% 5) Add 1.5 kg for sensors, hoses, and fittings

% Unpack Inputs
inj_material = IN.inj_material;
TCA_material = IN.TCA_material;


% Constants
flange_thickness = 0.5 * 0.0254; % m
regen_thickness = 0.275 * 0.025; % m
injector_height = 1.5 * 0.0254; % m
misc_mass = 1.5; % kg 

% Materials Dictionary
materials = ["Inconel", "SS316L", "Cu", "AlSi10Mg"];
densities = [8230, 8000, 8960, 2680]; % kg/m^3
material_dict = dictionary(materials, densities);

% TCA Mass
chamber_diameter = 3 * 0.0254; % m
csv_filepath = 'engine_contour.xlsx';  % Change this to your CSV file path
thickness = regen_thickness;     % Constant thickness

%% Read 2D line from CSV
fprintf('Reading 2D line from %s...\n', csv_filepath);
points_2d = readmatrix(csv_filepath);

% If there's a header, skip it
if isnan(points_2d(1,1))
    points_2d = points_2d(2:end, :);
end

% Take only first two columns (x, y)
points_2d = points_2d(:, 1:2);

fprintf('Loaded %d points\n', size(points_2d, 1));

% Volume calc with rectangular segments
fprintf('\nCalculating volume with thickness = %.4f...\n', thickness);

% Calculate distance between each consecutive point
segments = diff(points_2d);
segment_lengths = sqrt(sum(segments.^2, 2));

% Each segment is a rectangular prism with:
% - length = distance between points
% - width = thickness
% - height = thickness
% Volume of each segment = length * thickness^2
segment_volumes = segment_lengths * thickness^2;

% Total volume is sum of all segments
total_volume = sum(segment_volumes);
total_length = sum(segment_lengths);

TCA_mass = (2.5 + (total_volume * material_dict(TCA_material))) * 1.2;

% Injector Mass 
flange_R = chamber_diameter/2 + regen_thickness + flange_thickness;
inj_volume = (flange_R^2 * pi) * injector_height; % m^3
inj_mass = inj_volume * material_dict(inj_material);

% Actuators Mass
actuators_mass = 10; % kg
tvc_natural_freq = 1.8 / 0.01;

% Total Mass
engine_mass = inj_mass + TCA_mass + actuators_mass + misc_mass; % kg

%% =======================
% Pack Output
%% =======================
ENGINE = struct();
ENGINE.mass = engine_mass;
ENGINE.TCA_mass = TCA_mass;
ENGINE.TVC_mass = actuators_mass;
ENGINE.inj_mass = inj_mass;
ENGINE.tvc_natural_freq = tvc_natural_freq;

end
