%% Initial Parameters
tic;

% Initial state
x = 0; y = 0; z = 55;
vx = 5; vy = 0; vz = -10;

% Final state
xf = 50; yf = 0; zf = 0;
vxf = 0; vyf = 0; vzf = -0.1;

% Vehicle parameters
m0 = 113; % wet mass [kg]
Tmax = 2446; % maximum thrust [N]
specificImpulse = 225; % [s]

%% Setup for GFOLD
r0 = [x y z];
v0 = [vx vy vz];
rf = [xf yf zf];
vf = [vxf vyf vzf];

%% Call GFOLD solver
sol = gfold1(r0, v0, rf, vf, m0, Tmax, specificImpulse);

%% Display and plot results
if ~isempty(sol)
    disp('Solution found!');
    disp(sol);
    
    % Plot trajectory if solution exists
    plotTrajectory(sol);
else
    disp('No solution found');
end

totalTime = toc;
fprintf('Total computation time: %.2f seconds\n', totalTime);

%% Helper Functions
function sol = gfold1(r0, v0, rf, vf, m0, Tmax, Isp)
    % Set up Python environment
    pyenv('Version', 'C:\Users\joshu\AppData\Local\Programs\Python\Python311\python.exe');
    
    % Add the g-fold directory to Python path if needed
    % Only needed if not installed via pip
    % gfold_path = 'C:\path\to\g-fold\generator';
    % if count(py.sys.path, gfold_path) == 0
    %     insert(py.sys.path, int32(0), gfold_path);
    % end
    
    try
        % Import modules
        np = py.importlib.import_module('numpy');
        gfold = py.importlib.import_module('gfold');
        
        % Create solver (no parameters in constructor)
        solver = gfold.GFoldSolver();
        
        % Update parameters using the update_parameter method
        % Based on the code, these are the parameters that can be updated
        solver.update_parameter('m_wet', m0);
        solver.update_parameter('T_max', Tmax);
        solver.update_parameter('I_sp', Isp);
        
        % Set initial state
        solver.update_parameter('r_I_init', np.array(py.list(r0)));
        solver.update_parameter('v_I_init', np.array(py.list(v0)));
        
        % Set final state
        solver.update_parameter('r_I_final', np.array(py.list(rf)));
        solver.update_parameter('v_I_final', np.array(py.list(vf)));
        
        % Solve the problem
        result = solver.solve(pyargs('verbose', true));
        
        % Convert Python result to MATLAB structure
        sol = struct();
        
        % Extract results from the solution dictionary
        % The solution is a Python dict with keys like 'position', 'velocity', etc.
        if ~isempty(result)
            % Convert numpy arrays to MATLAB
            sol.position = double(result{'position'});  % r_I
            sol.velocity = double(result{'velocity'});  % v_I
            sol.thrust = double(result{'thrust'});      % T_B
            sol.mass = double(result{'mass'});          % m
            sol.time = double(result{'time'});          % t
            sol.final_mass = double(result{'final_mass'});
            sol.fuel_used = m0 - sol.final_mass;
            sol.status = 'optimal';
        else
            sol = [];
        end
        
    catch ME
        disp('Error calling Python GFOLD:');
        disp(ME.message);
        for i = 1:length(ME.stack)
            fprintf('  In %s at line %d\n', ME.stack(i).name, ME.stack(i).line);
        end
        sol = [];
    end
end

function plotTrajectory(sol)
    if isempty(sol)
        return;
    end
    
    figure('Position', [100 100 1400 900]);
    
    % 3D Trajectory
    subplot(2,3,1);
    plot3(sol.position(:,1), sol.position(:,2), sol.position(:,3), 'b-', 'LineWidth', 2);
    grid on;
    xlabel('X [m]'); ylabel('Y [m]'); zlabel('Z [m]');
    title('3D Trajectory');
    axis equal;
    hold on;
    plot3(sol.position(1,1), sol.position(1,2), sol.position(1,3), 'go', 'MarkerSize', 12, 'MarkerFaceColor', 'g');
    plot3(sol.position(end,1), sol.position(end,2), sol.position(end,3), 'ro', 'MarkerSize', 12, 'MarkerFaceColor', 'r');
    legend('Trajectory', 'Start', 'Landing', 'Location', 'best');
    view(3);
    
    % Position vs Time
    subplot(2,3,2);
    plot(sol.time, sol.position(:,1), 'r-', 'LineWidth', 1.5); hold on;
    plot(sol.time, sol.position(:,2), 'g-', 'LineWidth', 1.5);
    plot(sol.time, sol.position(:,3), 'b-', 'LineWidth', 1.5);
    grid on;
    xlabel('Time [s]'); ylabel('Position [m]');
    legend('X', 'Y', 'Z', 'Location', 'best');
    title('Position vs Time');
    
    % Velocity vs Time
    subplot(2,3,3);
    plot(sol.time, sol.velocity(:,1), 'r-', 'LineWidth', 1.5); hold on;
    plot(sol.time, sol.velocity(:,2), 'g-', 'LineWidth', 1.5);
    plot(sol.time, sol.velocity(:,3), 'b-', 'LineWidth', 1.5);
    grid on;
    xlabel('Time [s]'); ylabel('Velocity [m/s]');
    legend('V_x', 'V_y', 'V_z', 'Location', 'best');
    title('Velocity vs Time');
    
    % Thrust Magnitude
    subplot(2,3,4);
    thrust_mag = sqrt(sum(sol.thrust.^2, 2));
    plot(sol.time, thrust_mag, 'k-', 'LineWidth', 1.5);
    grid on;
    xlabel('Time [s]'); ylabel('Thrust [N]');
    title('Thrust Magnitude vs Time');
    yline(max(thrust_mag), 'r--', 'Max Thrust');
    
    % Mass vs Time
    subplot(2,3,5);
    plot(sol.time, sol.mass, 'm-', 'LineWidth', 1.5);
    grid on;
    xlabel('Time [s]'); ylabel('Mass [kg]');
    title('Mass vs Time');
    
    % Thrust Vector Components
    subplot(2,3,6);
    plot(sol.time, sol.thrust(:,1), 'r-', 'LineWidth', 1.5); hold on;
    plot(sol.time, sol.thrust(:,2), 'g-', 'LineWidth', 1.5);
    plot(sol.time, sol.thrust(:,3), 'b-', 'LineWidth', 1.5);
    grid on;
    xlabel('Time [s]'); ylabel('Thrust [N]');
    legend('T_x', 'T_y', 'T_z', 'Location', 'best');
    title('Thrust Components vs Time');
    
    % Print summary
    fprintf('\n=== GFOLD Solution Summary ===\n');
    fprintf('Flight time: %.2f s\n', sol.time(end));
    fprintf('Initial mass: %.2f kg\n', sol.mass(1));
    fprintf('Final mass: %.2f kg\n', sol.final_mass);
    fprintf('Fuel used: %.2f kg\n', sol.fuel_used);
    fprintf('Max thrust: %.2f N\n', max(thrust_mag));
    fprintf('Avg thrust: %.2f N\n', mean(thrust_mag));
    fprintf('Initial altitude: %.2f m\n', sol.position(1,3));
    fprintf('Final position: [%.2f, %.2f, %.2f] m\n', sol.position(end,1), sol.position(end,2), sol.position(end,3));
    fprintf('Final velocity: [%.2f, %.2f, %.2f] m/s\n', sol.velocity(end,1), sol.velocity(end,2), sol.velocity(end,3));
    fprintf('==============================\n\n');
end

%% 
%% Run this first to discover the correct parameter names

%% Run this first to discover the correct parameter names
%% Inspect the config object
pyenv('Version', 'C:\Users\joshu\AppData\Local\Programs\Python\Python311\python.exe');

try
    gfold = py.importlib.import_module('gfold');
    solver = gfold.GFoldSolver();
    
    % Look at the config object
    disp('=== Config Object ===');
    config = solver.config;
    disp(config);
    
    % Get config attributes
    disp('=== Config Attributes ===');
    config_dir = py.dir(config);
    config_list = string(config_dir);
    
    % Filter out private/magic methods (those starting with __)
    public_attrs = config_list(~startsWith(config_list, '__'));
    disp(public_attrs);
    
    % Try to print each attribute value
    disp('=== Config Values ===');
    for i = 1:length(public_attrs)
        attr = public_attrs(i);
        try
            value = py.getattr(config, attr);
            fprintf('%s = %s\n', attr, string(value));
        catch
            fprintf('%s = <not displayable>\n', attr);
        end
    end
    
    % Also check parameters
    disp('=== Parameters Object ===');
    params = solver.parameters;
    disp(params);
    
catch ME
    disp('Error:');
    disp(ME.message);
end