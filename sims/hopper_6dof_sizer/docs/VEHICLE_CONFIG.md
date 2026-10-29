# Vehicle config files

Each vehicle has one config file, the sim's source of truth for the vehicle itself: geometry, masses and hardware.

| Vehicle | File |
|---|---|
| hopper | `vehicles/hopper/hopper_config.m` |
| ctv | `vehicles/ctv/ctv_config.m` (starts from the hopper's file and overrides it) |

The `VEHICLE` variant axis picks the file ([VARIANTS.md](VARIANTS.md)). Revisions are tracked in git, so there is one file per vehicle.

What does **not** go in it: requirements, the reference trajectory, wind, Monte Carlo dispersions and FDIR settings. Those stay in `inputs/mission_inputs.m` (and `inputs/ctv_inputs.m` for the CTV's mission).

## Changing a value

Edit the config file and run the sim. Nothing else needs to change: every setup script and shared subsystem reads the value from there.

Two things are generated from the geometry and need regenerating by hand:
- **Baffle damping table.** After changing the tanks or baffles, run `generateBaffleLUT`. `sim_setup` warns when `data/baffle_lut.mat` was built for different tanks.
- **CG/MOI tables** are rebuilt every setup, so no action is needed.

## Frame

- **Body x** points up the thrust axis from the base reference point, the bottom centre of the vehicle (x = 0 on the pad).
- **y and z** complete a right-handed set.
- Every position in the file is an x in this frame, in metres.
- Legs and lidars are placed by azimuth about body x, measured from +y.

## Fields

`VC = vehicleConfig(vehicle)` reads the file, checks it, and adds the derived values. These checks raise an error naming the field:
- required fields are present
- sizes, masses, thrust and stiffness are positive
- the engine, gimbal, centre of pressure and tanks lie inside `height`
- the two tanks don't overlap. This is skipped for a vehicle that burns no propellant.

| Field | Units | Meaning |
|---|---|---|
| `name`, `revision`, `source` | text | which vehicle, which design point, where the numbers came from |
| `height` | m | base to top |
| `burns_propellant` | logical | default true; false for the CTV (constant mass) |
| `mass.dry` | kg | dry mass from the mass budget (legs and intertank included) |
| `mass.dry_cg`, `mass.dry_I` | m, kg·m² | optional CAD dry CG (axial) and inertia about it; without them the dry block comes from `components` |
| `tanks.radius`, `.geometry`, `.material` | m, text | inner radius, layout, material |
| `tanks.ox/fu.height`, `.cg_x`, `.wall` | m | as-built length, centre, wall thickness |
| `tanks.ox/fu.mass` | kg | empty tank |
| `tanks.ox/fu.volume`, `.prop_mass` | m³, kg | tank volume, propellant load |
| `tanks.ox/fu.density`, `.viscosity` | kg/m³, m²/s | propellant (slosh and baffle damping) |
| `tanks.ox/fu.baffle` | -, m, m | ring baffles: `rings`, `width`, `thickness` |
| `engine.mass`, `.cg_x`, `.gimbal_x` | kg, m, m | engine assembly, its CG, the gimbal pivot |
| `engine.offset` | m | thrust-line offset [y z] |
| `engine.nominal_thrust`, `.throttle_range` | N, fraction | thrust and its [min max] |
| `engine.throttle.*` | N/s, s, N | throttle valve `rate_limit`, `latency`, `resolution` |
| `engine.tvc.*` | rad, rad/s, s | `max_angle`, `max_rate`, `latency`, `rise_time`, `damping`, `resolution` |
| `rcs.*` | N, m, s, N·m | roll fans: `max_thrust`, moment `arm`, spin-up `tau`, `resolution`, `throttle_rate` |
| `legs.*` | -, m, deg, N/m, kg | `n`, `foot_radius`, `azimuth` (in leg order), `stiffness` per leg, `damping` (× √wet mass), `mass`, `intertank_mass` |
| `components.<name>` | kg, m | `mass`, `cg_x` of structures, avionics, fluids, payload and COPV, for the dry block until a CAD export |
| `structure.modal.*` | Hz, - | first bending mode: a stopgap until a structural response model exists |
| `aero.area`, `.cd`, `.cp_x` | m², -, m | side area, drag coefficient, centre of pressure |
| `sensors.lidar.*` | deg, m | mounting: `azimuth`, `cant`, `radius`, `mount_x` |

**Derived values (from `vehicleConfig`):**
- `tanks.ox/fu.floor_x`: the tank bottom, read by `Slosh` and `Mass_Properties`.
- `mass.wet`.
- `legs.r_foot`: 3 × n foot positions.
- `engine.tvc.natural_freq` = 1.8 / rise_time.
- `structure.modal.k_engine` = (engine cg_x / height)².

Don't set these in the file.

**Tags in the file:**
- **`sizer 2026-10`**: a value frozen from the conceptual sizer when it stopped running in sim setup, such as tank wall masses and leg stiffness. Replace it with an as-built value.
- **`ASSUMPTION`**: a placeholder that hasn't been measured.

## How the sim reads it

- **`mission_inputs`** calls `vehicleInputs(IN, vehicleConfig(hopperVehicle()))`, which:
  - copies the config into the `IN` fields the older code reads (`IN.engine_cg`, `IN.tanks.max_radius`, ...)
  - puts the config itself in `IN.vehicle`.
- **Monte Carlo dispersions** then change the `IN` fields. `IN.vehicle` always stays the nominal design.
- **`vehicleSetup(IN)`** builds the structs the models read: `TANKS`, `STRUCT`, `ENGINE`, `VEH`, `OUT`, `r_leg`. Both `sim_setup` (through `main.m`) and `mc_sim_setup` use it, so they always agree.
- **The shared subsystems** read `IN.vehicle` directly:
  - `Slosh` and `Mass_Properties`: the tank radius and floors
  - `EDF_Dynamics`: the fan arm
  - `Ground_Contact`: leg stiffness and feet.

New code should read `IN.vehicle.*`. Don't add more `IN` copies.

## The sizer

The sizer (`size_tanks`, `size_structures`, `size_engine`, `size_avionics`) is a design-trade tool now:

```matlab
sizeVehicle('hopper')                               % sizer vs config, side by side
sizeVehicle('hopper', 'Write', 'outputs/hopper_draft_config.m')   % draft config with the sized values
```

Review a draft and copy the values you adopt into the config file. `writeVehicleConfig` writes any config struct as a file at full precision, so a draft reads back exactly.

## CAD mass properties

`importCadMassProps` is a stub until the CAD mass-properties export exists. Its help text says what it will do, including the CAD-to-body frame conversion the CTV uses (`VC.cad.to_body`). Until then, type CAD values into `mass.dry`, `mass.dry_cg` and `mass.dry_I`.

## Adding a vehicle

1. Write `vehicles/<name>/<name>_config.m`. Start from the closest existing file, or call it and override fields, as `ctv_config` does.
2. Add `<name>` to `VEHICLE` in `src/setup/variantAxes.m`, and its variant defaults in `variantConfig.m`.
3. If its mission differs, add `inputs/<name>_inputs.m` and call it at the end of `mission_inputs`, as the CTV does.
4. Run `buildtool test`. `VehicleConfigTest` loads every listed config.
