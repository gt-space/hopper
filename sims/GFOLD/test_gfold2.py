import sys
sys.path.append(r"C:\Users\joshu")

from gfold import GFoldSolver
import numpy as np
import matplotlib.pyplot as plt

# Create solver and solve
solver = GFoldSolver()
result = solver.solve(verbose=False)

# Print key results
print("=== GFOLD Solution Summary ===")
print(f"Optimal fuel cost: {result['solution_value']:.4f}")
print(f"Final mass: {result['final_mass']:.4f} kg")
print(f"Number of time points: {len(result['time_points'])}")
print(f"Mission duration: {result['time_points'][-1]:.2f} seconds")

# Check final state
final_pos = result['positions'][-1, :]  # Fixed indexing
final_vel = result['velocities'][-1, :]
print(f"\nFinal position: {final_pos}")
print(f"Final velocity: {final_vel}")

# Trajectory shape
print(f"\nTrajectory shape:")
print(f"  Positions: {result['positions'].shape}")
print(f"  Velocities: {result['velocities'].shape}")
print(f"  Thrusts: {result['thrusts'].shape}")

# Plot the trajectory
fig, axes = plt.subplots(2, 2, figsize=(12, 10))

# Position plot (X-Z plane)
axes[0, 0].plot(result['positions'][:, 0], result['positions'][:, 2])
axes[0, 0].scatter(result['positions'][0, 0], result['positions'][0, 2], 
                   c='green', s=100, marker='o', label='Start', zorder=5)
axes[0, 0].scatter(result['positions'][-1, 0], result['positions'][-1, 2], 
                   c='red', s=100, marker='x', label='Landing', zorder=5)
axes[0, 0].set_xlabel('X position (m)')
axes[0, 0].set_ylabel('Z position (m)')
axes[0, 0].set_title('Trajectory (X-Z plane)')
axes[0, 0].legend()
axes[0, 0].grid(True)

# Velocity plot
axes[0, 1].plot(result['time_points'], result['velocities'][:, 0], label='Vx')
axes[0, 1].plot(result['time_points'], result['velocities'][:, 1], label='Vy')
axes[0, 1].plot(result['time_points'], result['velocities'][:, 2], label='Vz')
axes[0, 1].set_xlabel('Time (s)')
axes[0, 1].set_ylabel('Velocity (m/s)')
axes[0, 1].set_title('Velocity vs Time')
axes[0, 1].legend()
axes[0, 1].grid(True)

# Thrust magnitude
thrust_mag = result['thrusts']  # Already magnitude
axes[1, 0].plot(result['time_points'], thrust_mag)
axes[1, 0].set_xlabel('Time (s)')
axes[1, 0].set_ylabel('Thrust (N)')
axes[1, 0].set_title('Thrust Magnitude')
axes[1, 0].grid(True)

# 3D trajectory
from mpl_toolkits.mplot3d import Axes3D
axes[1, 1].remove()
ax3d = fig.add_subplot(2, 2, 4, projection='3d')
ax3d.plot(result['positions'][:, 0], result['positions'][:, 1], result['positions'][:, 2])
ax3d.scatter(result['positions'][0, 0], result['positions'][0, 1], result['positions'][0, 2],
             c='green', s=100, marker='o', label='Start')
ax3d.scatter(result['positions'][-1, 0], result['positions'][-1, 1], result['positions'][-1, 2],
             c='red', s=100, marker='x', label='Landing')
ax3d.set_xlabel('X (m)')
ax3d.set_ylabel('Y (m)')
ax3d.set_zlabel('Z (m)')
ax3d.set_title('3D Trajectory')
ax3d.legend()

plt.tight_layout()
plt.savefig('gfold_trajectory.png', dpi=150)
plt.show()

print("\nPlot saved as 'gfold_trajectory.png'")

# Print initial and final conditions
print("\n=== Initial Conditions ===")
print(f"Position: {result['positions'][0, :]}")
print(f"Velocity: {result['velocities'][0, :]}")

print("\n=== Final Conditions ===")
print(f"Position: {result['positions'][-1, :]}")
print(f"Velocity: {result['velocities'][-1, :]}")