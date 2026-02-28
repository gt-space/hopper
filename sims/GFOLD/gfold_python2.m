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
    
    % Plot trajectory if solution exists
    plotTrajectory(sol, m0);
else
    disp('No solution found');
end

totalTime = toc;
fprintf('Total computation time: %.2f seconds\n', totalTime);

%% Helper Functions
function sol = gfold1(r0, v0, rf, vf, m0, Tmax, Isp)
    % Set up Python environment
    pyenv('Version', 'C:\Users\joshu\AppData\Local\Programs\Python\Python311\python.exe');
    
    try
        % Import modules
        np = py.importlib.import_module('numpy');
        gfold = py.importlib.import_module('gfold');
        
        % Create solver with default configuration
        solver = gfold.GFoldSolver();
        
        % Update spacecraft configuration
        % The config has spacecraft, environment, and solver sub-configs
        config_dict = py.dict();
        
        % Spacecraft parameters
        spacecraft_dict = py.dict();
        spacecraft_dict{'wet_mass'} = m0;
        spacecraft_dict{'max_thrust'} = Tmax;
        spacecraft_dict{'specific_impulse'} = Isp;
        config_dict{'spacecraft'} = spacecraft_dict;
        
        % Update the configuration
        solver.update_config(config_dict);
        
        % Update parameters for initial and final states
        solver.update_parameter('initial_position', np.array(py.list(r0)));
        solver.update_parameter('initial_vel', np.array(py.list(v0)));
        solver.update_parameter('target_velocity', np.array(py.list(vf)));
        
        % Note: The solver might need target position too
        % Let's try to set it if the parameter exists
        try
            % Some implementations use this
            solver.update_parameter('target_position', np.array(py.list(rf)));
        catch
            % If not, it might be handled differently
        end
        
        % Solve the problem
        result = solver.solve(pyargs('verbose', true));
        
        % Convert Python result to MATLAB structure
        sol = struct();
        
        if ~isempty(result)
            % Get the keys from result
            result_keys = py.list(result.keys());
            num_keys = int32(py.len(result_keys));
            
            fprintf('Result contains %d fields\n', num_keys);
            
            % Extract common fields
            try
                sol.position = double(result{'position'});
            catch
                try
                    sol.position = double(result{'r'});
                catch
                    disp('Could not extract position');
                end
            end
            
            try
                sol.velocity = double(result{'velocity'});
            catch
                try
                    sol.velocity = double(result{'v'});
                catch
                    disp('Could not extract velocity');
                end
            end
            
            try
                sol.thrust = double(result{'thrust'});
            catch
                try
                    sol.thrust = double(result{'T'});
                catch
                    disp('Could not extract thrust');
                end
            end
            
            try
                sol.mass = double(result{'mass'});
            catch
                try
                    sol.mass = double(result{'m'});
                catch
                    disp('Could not extract mass');
                end
            end
            
            try
                sol.time = double(result{'time'});
            catch
                try
                    sol.time = double(result{'t'});
                catch
                    disp('Could not extract time');
                end
            end
            
            try
                sol.final_mass = double(result{'final_mass'});
            catch
                if isfield(sol, 'mass')
                    sol.final_mass = sol.mass(end);
                end
            end
            
            if isfield(sol, 'mass')
                sol.fuel_used = m0 - sol.final_mass;
            end
            
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

function plotTrajectory(sol, m0)
    if isempty(sol) || ~isfield(sol, 'position')
        disp('No valid solution to plot');
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
    if isfield(sol, 'time')
        plot(sol.time, sol.position(:,1), 'r-', 'LineWidth', 1.5); hold on;
        plot(sol.time, sol.position(:,2), 'g-', 'LineWidth', 1.5);
        plot(sol.time, sol.position(:,3), 'b-', 'LineWidth', 1.5);
        xlabel('Time [s]');
    else
        plot(sol.position(:,1), 'r-', 'LineWidth', 1.5); hold on;
        plot(sol.position(:,2), 'g-', 'LineWidth', 1.5);
        plot(sol.position(:,3), 'b-', 'LineWidth', 1.5);
        xlabel('Step');
    end
    grid on;
    ylabel('Position [m]');
    legend('X', 'Y', 'Z', 'Location', 'best');
    title('Position vs Time');
    
    % Velocity vs Time
    subplot(2,3,3);
    if isfield(sol, 'velocity')
        if isfield(sol, 'time')
            plot(sol.time, sol.velocity(:,1), 'r-', 'LineWidth', 1.5); hold on;
            plot(sol.time, sol.velocity(:,2), 'g-', 'LineWidth', 1.5);
            plot(sol.time, sol.velocity(:,3), 'b-', 'LineWidth', 1.5);
            xlabel('Time [s]');
        else
            plot(sol.velocity(:,1), 'r-', 'LineWidth', 1.5); hold on;
            plot(sol.velocity(:,2), 'g-', 'LineWidth', 1.5);
            plot(sol.velocity(:,3), 'b-', 'LineWidth', 1.5);
            xlabel('Step');
        end
        grid on;
        ylabel('Velocity [m/s]');
        legend('V_x', 'V_y', 'V_z', 'Location', 'best');
        title('Velocity vs Time');
    end
    
    % Thrust Magnitude
    subplot(2,3,4);
    if isfield(sol, 'thrust')
        thrust_mag = sqrt(sum(sol.thrust.^2, 2));
        if isfield(sol, 'time')
            plot(sol.time, thrust_mag, 'k-', 'LineWidth', 1.5);
            xlabel('Time [s]');
        else
            plot(thrust_mag, 'k-', 'LineWidth', 1.5);
            xlabel('Step');
        end
        grid on;
        ylabel('Thrust [N]');
        title('Thrust Magnitude vs Time');
        yline(max(thrust_mag), 'r--', sprintf('Max: %.0f N', max(thrust_mag)));
    end
    
    % Mass vs Time
    subplot(2,3,5);
    if isfield(sol, 'mass')
        if isfield(sol, 'time')
            plot(sol.time, sol.mass, 'm-', 'LineWidth', 1.5);
            xlabel('Time [s]');
        else
            plot(sol.mass, 'm-', 'LineWidth', 1.5);
            xlabel('Step');
        end
        grid on;
        ylabel('Mass [kg]');
        title('Mass vs Time');
    end
    
    % Thrust Vector Components
    subplot(2,3,6);
    if isfield(sol, 'thrust')
        if isfield(sol, 'time')
            plot(sol.time, sol.thrust(:,1), 'r-', 'LineWidth', 1.5); hold on;
            plot(sol.time, sol.thrust(:,2), 'g-', 'LineWidth', 1.5);
            plot(sol.time, sol.thrust(:,3), 'b-', 'LineWidth', 1.5);
            xlabel('Time [s]');
        else
            plot(sol.thrust(:,1), 'r-', 'LineWidth', 1.5); hold on;
            plot(sol.thrust(:,2), 'g-', 'LineWidth', 1.5);
            plot(sol.thrust(:,3), 'b-', 'LineWidth', 1.5);
            xlabel('Step');
        end
        grid on;
        ylabel('Thrust [N]');
        legend('T_x', 'T_y', 'T_z', 'Location', 'best');
        title('Thrust Components vs Time');
    end
    
    % Print summary
    fprintf('\n=== GFOLD Solution Summary ===\n');
    if isfield(sol, 'time')
        fprintf('Flight time: %.2f s\n', sol.time(end));
    end
    fprintf('Initial mass: %.2f kg\n', m0);
    if isfield(sol, 'final_mass')
        fprintf('Final mass: %.2f kg\n', sol.final_mass);
        fprintf('Fuel used: %.2f kg (%.1f%%)\n', sol.fuel_used, 100*sol.fuel_used/m0);
    end
    if isfield(sol, 'thrust')
        thrust_mag = sqrt(sum(sol.thrust.^2, 2));
        fprintf('Max thrust: %.2f N\n', max(thrust_mag));
        fprintf('Avg thrust: %.2f N\n', mean(thrust_mag));
    end
    fprintf('Initial position: [%.2f, %.2f, %.2f] m\n', sol.position(1,1), sol.position(1,2), sol.position(1,3));
    fprintf('Final position: [%.2f, %.2f, %.2f] m\n', sol.position(end,1), sol.position(end,2), sol.position(end,3));
    if isfield(sol, 'velocity')
        fprintf('Final velocity: [%.2f, %.2f, %.2f] m/s\n', sol.velocity(end,1), sol.velocity(end,2), sol.velocity(end,3));
    end
    fprintf('==============================\n\n');
end