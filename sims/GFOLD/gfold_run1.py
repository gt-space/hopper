import sys
sys.path.append(r"C:\Users\joshu")

from gfold import GFoldSolver
import numpy as np
import matplotlib.pyplot as plt

# ====================================
# YOUR VEHICLE PARAMETERS - ASCENT-HOVER-DESCENT MISSION
# ====================================

# IMPORTANT NOTE: GFOLD is designed for powered descent (landing) optimization.
# It minimizes fuel while landing from an initial state to a target.
# For your ascent-hover-descent mission, we'll need to break it into phases:
# 
# Option 1: Run GFOLD for descent phase only (from 55m hover to landing)
# Option 2: Use GFOLD creatively by setting "initial" as top of trajectory

# Let's start with DESCENT PHASE (from 55m hover to ground)
# This is what GFOLD is designed for!

custom_config = {
    # Spacecraft mass parameters
    'wet_mass': 113.75,            # Total mass with full fuel (kg)
    'fuel': 18.83,                 # Total propellant: LOX 16kg + Dodecane 2.83kg
                                   # NOTE: If total is 25kg (16+9), change to 'fuel': 25.0
    
    # Thrust parameters
    'real_max_thrust': 2446.521,   # Maximum thrust (N)
    'min_thrust_pct': 0.3636,      # 889.644 / 2446.521
    'max_thrust_pct': 1.0,         # Using full capability
    
    # Fuel consumption rate
    # Estimate from Isp: For LOX/Dodecane Isp ~280-300s
    # fuel_consumption = 1/(Isp * g0) = 1/(290 * 9.81) ≈ 0.000352
    'fuel_consumption': 0.00035,   # kg/N/s
    
    # DESCENT PHASE: Starting from hover at 55m
    'initial_position': [0, 0, 55],        # Hovering at 55m altitude
    'initial_velocity': [0, 0, 0],         # Stationary hover (zero velocity)
    
    # Target state (landing)
    'target_position': [0, 0, 0],          # Landing at ground level
    'target_velocity': [0, 0, -1.25],      # Soft landing: 1.25 m/s downward (~4 ft/s)
    
    # Environment - EARTH GRAVITY
    'gravity': [0, 0, -9.81],       # Earth gravity
    
    # Constraints
    'glide_slope_angle': 0,         # No glide slope constraint
    'max_angle': 20,                # Keep thrust mostly vertical (20 deg from vertical)
    'max_velocity': 50,             # Maximum velocity constraint (m/s)
    
    # Solver settings
    'n': 80,                        # Number of discretization points
    'time_of_flight': 15.0,         # Estimated descent time from 55m hover (seconds)
}

# ====================================
# VERIFY CONFIGURATION
# ====================================
dry_mass = custom_config['wet_mass'] - custom_config['fuel']
print("=== Vehicle Configuration ===")
print(f"Dry Mass: {dry_mass:.2f} kg")
print(f"Wet Mass: {custom_config['wet_mass']:.2f} kg")
print(f"Propellant Mass: {custom_config['fuel']:.2f} kg ({100*custom_config['fuel']/custom_config['wet_mass']:.2f}%)")
print(f"Min Thrust: {custom_config['real_max_thrust'] * custom_config['min_thrust_pct']:.3f} N")
print(f"Max Thrust: {custom_config['real_max_thrust']:.3f} N")
print(f"Initial TWR (wet): {custom_config['real_max_thrust'] / (custom_config['wet_mass'] * 9.81):.2f}")
print(f"Final TWR (dry): {custom_config['real_max_thrust'] / (dry_mass * 9.81):.2f}")
print(f"\n=== Mission Profile ===")
print(f"Phase: DESCENT from hover")
print(f"Starting altitude: {custom_config['initial_position'][2]} m")
print(f"Target landing velocity: {abs(custom_config['target_velocity'][2]):.2f} m/s")
print()

# ====================================
# ESTIMATE FUEL FOR HOVER (MANUAL CALCULATION)
# ====================================
# For hovering, thrust must equal weight
# After ascent, let's estimate mass = 110 kg (used ~3.75 kg fuel in ascent)
hover_mass = 110.0  # kg (estimate after ascent)
hover_thrust = hover_mass * 9.81  # N
hover_duration = 3.0  # seconds

# Fuel consumption = fuel_consumption * thrust * time
hover_fuel = custom_config['fuel_consumption'] * hover_thrust * hover_duration
print(f"=== Hover Phase Estimate (manual) ===")
print(f"Hover mass: {hover_mass:.2f} kg")
print(f"Hover thrust required: {hover_thrust:.2f} N")
print(f"Hover duration: {hover_duration:.1f} s")
print(f"Estimated hover fuel: {hover_fuel:.2f} kg")
print()

# Adjust wet mass for descent (after ascent + hover fuel burned)
estimated_ascent_fuel = 3.75  # kg (rough estimate)
descent_wet_mass = custom_config['wet_mass'] - estimated_ascent_fuel - hover_fuel
print(f"Estimated mass at start of descent: {descent_wet_mass:.2f} kg")

# Update config for descent phase
custom_config['wet_mass'] = descent_wet_mass
custom_config['fuel'] = custom_config['fuel'] - estimated_ascent_fuel - hover_fuel

print(f"Fuel available for descent: {custom_config['fuel']:.2f} kg")
print()

# ====================================
# SOLVE DESCENT PHASE
# ====================================

# Create solver with custom config
solver = GFoldSolver(config=custom_config)

# Solve the problem
print("=== Solving GFOLD for DESCENT PHASE ===")
result = solver.solve(verbose=False)

# ====================================
# PRINT RESULTS
# ====================================

print("\n=== GFOLD Descent Solution ===")
print(f"Optimal fuel cost: {result['solution_value']:.4f}")
print(f"Final mass: {result['final_mass']:.4f} kg")
print(f"Descent fuel used: {custom_config['wet_mass'] - result['final_mass']:.4f} kg")
print(f"Descent duration: {result['time_points'][-1]:.2f} seconds")

print(f"\n=== Initial Conditions (Start of Descent) ===")
print(f"Position: {result['positions'][0, :]} m")
print(f"Velocity: {result['velocities'][0, :]} m/s")

print(f"\n=== Final Conditions (Landing) ===")
print(f"Position: {result['positions'][-1, :]} m")
print(f"Velocity: {result['velocities'][-1, :]} m/s")
landing_vel = np.linalg.norm(result['velocities'][-1, :])
print(f"Landing velocity magnitude: {landing_vel:.3f} m/s")

# Check landing velocity
if landing_vel <= 1.35:  # Allow small margin
    print(f"✓ Landing velocity OK: {landing_vel:.3f} m/s ≈ 1.25 m/s target")
else:
    print(f"✗ Landing velocity: {landing_vel:.3f} m/s (target: 1.25 m/s)")

# Total mission fuel estimate
total_fuel_used = estimated_ascent_fuel + hover_fuel + (custom_config['wet_mass'] - result['final_mass'])
print(f"\n=== TOTAL MISSION FUEL ESTIMATE ===")
print(f"Ascent fuel (estimate): {estimated_ascent_fuel:.2f} kg")
print(f"Hover fuel (estimate): {hover_fuel:.2f} kg")
print(f"Descent fuel (GFOLD): {custom_config['wet_mass'] - result['final_mass']:.2f} kg")
print(f"TOTAL: {total_fuel_used:.2f} kg")
print(f"Fuel margin: {18.83 - total_fuel_used:.2f} kg ({100*(18.83 - total_fuel_used)/18.83:.1f}%)")

# ====================================
# PLOT RESULTS
# ====================================

fig, axes = plt.subplots(2, 3, figsize=(16, 10))

# Altitude vs Time (DESCENT ONLY)
axes[0, 0].plot(result['time_points'], result['positions'][:, 2], 'b-', linewidth=2.5)
axes[0, 0].scatter(0, 55, c='green', s=150, marker='o', label='Start Descent', 
                   zorder=5, edgecolors='black')
axes[0, 0].scatter(result['time_points'][-1], 0, c='red', s=150, marker='x', 
                   label='Landing', zorder=5, linewidths=3)
axes[0, 0].set_xlabel('Time (s)', fontsize=11)
axes[0, 0].set_ylabel('Altitude (m)', fontsize=11)
axes[0, 0].set_title('Descent Altitude vs Time', fontsize=12, fontweight='bold')
axes[0, 0].legend()
axes[0, 0].grid(True, alpha=0.3)
axes[0, 0].set_ylim(bottom=-2)

# Velocity plot
axes[0, 1].plot(result['time_points'], result['velocities'][:, 0], label='Vx', linewidth=2)
axes[0, 1].plot(result['time_points'], result['velocities'][:, 1], label='Vy', linewidth=2)
axes[0, 1].plot(result['time_points'], result['velocities'][:, 2], label='Vz (vertical)', linewidth=2)
axes[0, 1].axhline(y=-1.25, color='red', linestyle='--', alpha=0.5, label='Target Vz')
axes[0, 1].set_xlabel('Time (s)', fontsize=11)
axes[0, 1].set_ylabel('Velocity (m/s)', fontsize=11)
axes[0, 1].set_title('Velocity vs Time', fontsize=12, fontweight='bold')
axes[0, 1].legend()
axes[0, 1].grid(True, alpha=0.3)

# Thrust magnitude
thrust_mag = result['thrusts']
axes[0, 2].plot(result['time_points'], thrust_mag, 'b-', linewidth=2.5)
axes[0, 2].axhline(y=2446.521, color='r', linestyle='--', 
                   label=f'Max Thrust (2446.5 N)', alpha=0.7)
axes[0, 2].axhline(y=889.644, color='orange', linestyle='--', 
                   label=f'Min Thrust (889.6 N)', alpha=0.7)
# Add hover thrust reference
axes[0, 2].axhline(y=hover_thrust, color='green', linestyle=':', 
                   label=f'Hover Thrust (~{hover_thrust:.0f} N)', alpha=0.7)
axes[0, 2].set_xlabel('Time (s)', fontsize=11)
axes[0, 2].set_ylabel('Thrust (N)', fontsize=11)
axes[0, 2].set_title('Thrust Magnitude (Descent)', fontsize=12, fontweight='bold')
axes[0, 2].legend(fontsize=9)
axes[0, 2].grid(True, alpha=0.3)

# Vertical velocity vs altitude
axes[1, 0].plot(result['velocities'][:, 2], result['positions'][:, 2], 'b-', linewidth=2.5)
axes[1, 0].scatter(0, 55, c='green', s=150, marker='o', label='Start', 
                   zorder=5, edgecolors='black')
axes[1, 0].scatter(result['velocities'][-1, 2], 0, c='red', s=150, marker='x', 
                   label='Landing', zorder=5, linewidths=3)
axes[1, 0].axvline(x=-1.25, color='red', linestyle='--', alpha=0.5, label='Target landing Vz')
axes[1, 0].set_xlabel('Vertical Velocity (m/s)', fontsize=11)
axes[1, 0].set_ylabel('Altitude (m)', fontsize=11)
axes[1, 0].set_title('Descent Profile', fontsize=12, fontweight='bold')
axes[1, 0].legend()
axes[1, 0].grid(True, alpha=0.3)

# Mass over time
mass_over_time = custom_config['wet_mass'] - (custom_config['wet_mass'] - result['final_mass']) * (result['time_points'] / result['time_points'][-1])
axes[1, 1].plot(result['time_points'], mass_over_time, 'g-', linewidth=2.5)
axes[1, 1].axhline(y=dry_mass, color='r', linestyle='--', 
                   label=f'Dry Mass ({dry_mass:.1f} kg)', alpha=0.7)
axes[1, 1].set_xlabel('Time (s)', fontsize=11)
axes[1, 1].set_ylabel('Mass (kg)', fontsize=11)
axes[1, 1].set_title('Vehicle Mass vs Time (Descent)', fontsize=12, fontweight='bold')
axes[1, 1].legend()
axes[1, 1].grid(True, alpha=0.3)

# Trajectory (X-Z plane)
axes[1, 2].plot(result['positions'][:, 0], result['positions'][:, 2], 'b-', linewidth=2.5)
axes[1, 2].scatter(result['positions'][0, 0], result['positions'][0, 2], 
                   c='green', s=150, marker='o', label='Start', zorder=5, edgecolors='black')
axes[1, 2].scatter(result['positions'][-1, 0], result['positions'][-1, 2], 
                   c='red', s=150, marker='x', label='Landing', zorder=5, linewidths=3)
axes[1, 2].set_xlabel('X position (m)', fontsize=11)
axes[1, 2].set_ylabel('Z position (m)', fontsize=11)
axes[1, 2].set_title('Trajectory (X-Z plane)', fontsize=12, fontweight='bold')
axes[1, 2].legend()
axes[1, 2].grid(True, alpha=0.3)
axes[1, 2].set_ylim(bottom=-2)

plt.tight_layout()
plt.savefig('gfold_descent_trajectory.png', dpi=150, bbox_inches='tight')
plt.show()

print("\nPlot saved as 'gfold_descent_trajectory.png'")

# ====================================
# SAVE TRAJECTORY DATA
# ====================================
np.savez('descent_trajectory.npz',
         time=result['time_points'],
         position=result['positions'],
         velocity=result['velocities'],
         thrust=result['thrusts'],
         mass=mass_over_time)
print("Trajectory data saved as 'descent_trajectory.npz'")