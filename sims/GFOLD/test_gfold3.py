import sys
sys.path.append(r"C:\Users\joshu")

from gfold import GFoldSolver

# Create solver
solver = GFoldSolver()

print("=== Current GFoldSolver Configuration ===\n")

# Show all config attributes and their current values
print("Configuration attributes:")
config_attrs = ['wet_mass', 'log_dry_mass', 'log_wet_mass', 'max_thrust', 
                'min_thrust', 'real_max_thrust', 'cos_max_angle', 
                'sin_glide_slope', 'time_of_flight']

for attr in config_attrs:
    if hasattr(solver.config, attr):
        value = getattr(solver.config, attr)
        print(f"  {attr}: {value}")

# Check spacecraft sub-config
print("\nSpacecraft configuration:")
if hasattr(solver.config, 'spacecraft'):
    print(f"  Type: {type(solver.config.spacecraft)}")
    for attr in dir(solver.config.spacecraft):
        if not attr.startswith('_'):
            value = getattr(solver.config.spacecraft, attr)
            print(f"    {attr}: {value}")

# Check environment sub-config
print("\nEnvironment configuration:")
if hasattr(solver.config, 'environment'):
    print(f"  Type: {type(solver.config.environment)}")
    for attr in dir(solver.config.environment):
        if not attr.startswith('_'):
            value = getattr(solver.config.environment, attr)
            print(f"    {attr}: {value}")

# Check solver sub-config
print("\nSolver configuration:")
if hasattr(solver.config, 'solver'):
    print(f"  Type: {type(solver.config.solver)}")
    for attr in dir(solver.config.solver):
        if not attr.startswith('_'):
            value = getattr(solver.config.solver, attr)
            print(f"    {attr}: {value}")

# Check the to_dict method
print("\n=== Full Config Dictionary ===")
print(solver.config.to_dict())

# Check update methods
print("\n=== Available Update Methods ===")
print("update_config signature:", inspect.signature(solver.update_config) if 'inspect' in dir() else "import inspect first")
print("update_parameter signature:", inspect.signature(solver.update_parameter) if 'inspect' in dir() else "import inspect first")