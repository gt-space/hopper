function AVI = size_avionics(IN)
%% Description
% Power is drawn from solenoid valves, 2X TVC Actuators, 2X Throttle
% Valves, and avionics boards.

% Standby and Flight Time are separate hours cases that power must be
% supplied for when umbilical is dicsonnected and it depends on each component.

% Empirical correlations will be used to return the 

% Capabilities
% 1) Calculate required battery, boards, and harness mass
% 2) Calculate total power draw for the vehicle

% Battery Specs from here - https://www.18650batterystore.com/collections/18650-batteries

%% Unpack Inputs
standby_hours = IN.avionics.standby_hours; % hr
flight_time = IN.avionics.flight_time / 3600; % sec -> hr
lcoms_time = IN.avionics.lcoms_time / 3600; % sec -- > hr
battery_voltage = IN.avionics.voltage; % V

% Power Calculations
TVC_power = IN.tvc.voltage * IN.tvc.current; % W
Motor_power = IN.throttle_valve.voltage * IN.throttle_valve.current; % W
Boards_power = IN.avionics.boards.voltage * IN.avionics.boards.current; % W
Valves_power = IN.solenoid_valve_quantity * IN.valves.voltage * IN.valves.current; % W

% Capacity Calculations
energy_standby = standby_hours * (Boards_power + Motor_power); % Wh
energy_flight = flight_time * (TVC_power + Motor_power + Valves_power + Boards_power); % Wh
energy_lcoms = lcoms_time * Valves_power; % Wh
tot_load = energy_standby + energy_flight + energy_lcoms; % Wh
capacity = (tot_energy * (1 * 1 * 1.15)) / (48 * 0.2); % Ah

% Cell Sizing
V_load_max = 30 / 100;
V_charging = 4.2; % V
num_cells = ceil((battery_voltage * (1 + V_load_max)) / V_charging);

% Battery Specs (from website)
m_batt = 50 / 1000; % g --> kg

% Mass Tracker
boards_mass = 8 * 0.453592; % lbm --> kg, includes structural mass, 11 lbm for Vespula as baseline
battery_mass = m_batt * num_cells; % kg
harnesses_mass = 3 * 0.453592; % lbm --> kg
avi_mass = boards_mass + battery_mass + harnesses_mass; % kg

% Pack Outputs
AVI = struct();
AVI.tot_load = tot_load;
AVI.battery_capacity = capacity;
AVI.num_cells = num_cells;
AVI.boards_mass = boards_mass;
AVI.battery_mass = battery_mass;
AVI.harnesses_mass = harnesses_mass;
AVI.mass = avi_mass;

end