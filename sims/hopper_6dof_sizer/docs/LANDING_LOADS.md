# Landing loads model

`src/analysis/landing` computes the loads the legs see at touchdown and sizes crush cores. It is a standalone model, separate from the flight sim: it starts at first contact from a touchdown state and runs through settling or tip-over. The flight model's `Ground_Contact` still uses the sizer's simple feet (see *Not yet done*).

```matlab
G = landingGear(Core=crushCore('pad', Plateau=5e3));   % as-drawn legs + a crush pad
S = landingSurface('concrete');                         % or soil_firm, soil_soft, rigid
C = landingCase(SinkRate=3.048, Lateral=0.61, Heading='leg', Tilt=5);
R = landingSim(G, S, C);   % R.peak: strut loads, stroke, g, tilt; R.hist: time history

T  = landingSweep();       % cores x surfaces x cases -> results/landing_sweep.mat
Sm = landingReport(T);     % margins per core -> outputs/landing/*.csv, landing_sweep.png
```

### Crush-pad × engine-cutoff trade

```matlab
R = runMonteCarlo(200, UseParallel=true, OutputFile=hopperRoot('results','mc_touchdown_200.mat'));
P = load(hopperRoot('results','fluid_shutdown_profile.mat'));   % thrust after the cutoff command, fluid network
X = cutoffTrade(hopperRoot('results','mc_touchdown_200.mat'), ShutdownProfile=P.prof, Plateau=[4 5]*1e3);
```

- `touchdownCases` turns each flight run's last state into a `landingCase` (drift, tilt, rates, mass, thrust).
- `cutoffTrade` finds, per run, pad and surface, the highest cut height that doesn't crush the pads (or, without a profile, the crush-onset sink rate), and whether an engine-on landing crushes. Then it samples the cut-height error σ to get P(no crush) per commanded height. Each search is appended to a CSV log as it finishes, and a rerun skips what's logged.
- The shutdown curve came from flying `PROPULSION=fluid_network` with the cutoff at 1 m: measured to 0.5 s after the command, the tail extrapolated.

## What is modelled

- **Legs, as drawn.** Strut tube ends come from `Updated Vehicle Assy.STL` (CAD 2026-10-29) and are stored in `VC.legs.geom` in `hopper_config.m`; tube sections were measured from the same file:
  - main: 1.5/1.0 in OD/ID
  - aux: 1.0/0.75 in OD/ID
  - Each foot node (pad, fitting, half the strut mass) is joined to the body by three axial struts with exact lengths.
  - `Joints='fixed'` adds strut bending (clamped struts, ball-jointed pad) as a stiff bound.
- **Crush cores.** Elastic up to the peak load, then plateau, then densification at 75 % of the core length; the crush is permanent.
  - `'main'`: inline in each main strut.
  - `'pad'`: between the foot and the ground.
  - Aux struts run in tension on landing, so a compression core there never crushes.
- **Ground.**
  - Linear spring-damper for `rigid` and `concrete`.
  - Elastic-plastic subgrade for `soil_firm` and `soil_soft`: loading is ks × pad area, unloading is 10× stiffer, and the pad leaves a dent.
  - Friction is a stick-slip spring capped at μ × N.
- **Touchdown cases** (`landingCases`): flat, plus sideways drift toward a leg (Ansys "1-2-1": one foot, then two, then the last) or between two legs ("2-2"). Each drift case runs in two attitudes:
  - `'braking'`: the top is tilted away from the drift, as when killing sideways speed, so the trailing foot lands first and the vehicle rocks forward onto the leading foot;
  - `'leading'`: the top is tilted into the drift, so the leading foot lands first and digs in.
- **Body.** A rigid 6-DOF body. Mass, CG and inertia default to the 145 kg design mass of the leg hand calcs, with the CAD dry CG and inertia scaled to that mass. Thrust at contact defaults to 0, since the engine cuts off 5 cm up.

## Validation

`tests/unit/LandingModelTest.m`:
- **Nominal drop vs Ansys:** 145 kg at 1.5 m/s on rigid, frictionless ground gives main 42.8 kN and aux −18.8 kN. Structures' Ansys Flat11 case gives 45.2 / −18.1 kN.
- **Geometry vs CAD:** strut lengths and the static force ratios match the CAD.
- **Crush caps:** core and pad forces never exceed their plateau.
- **Energy:** no run ends with more energy than it started with.
- **Step size:** halving dt changes peaks by under 3 %.

## Things the model shows about the leg

- **The struts of a leg lie within 10° of each other.** As a pinned truss the leg is soft across the strut bundle (about 0.9 MN/m vertical), so:
  - the main strut carries about 5× the vertical foot load in compression;
  - the aux struts carry about 2.2× in tension.
- **Main-strut buckling, with the correct area moment π/64 (D⁴ − d⁴), is Pcr ≈ 69 kN** at the 35.7 in tube-end length. The leg hand calcs use π/32, the polar moment, which doubles Pcr.
- **Aux-strut buckling is Pcr ≈ 11.5 kN.** Lateral drift puts the aux struts in compression, and without a core they reach that load in off-nominal cases.

## Assumptions to replace

- Pad and fitting mass: 0.5 kg per foot (`VC.legs.foot_fitting_mass`).
- Surface stiffness, damping and friction are textbook values, not site data.
- The ground force acts at the foot node; the pad pivots.
- Joint friction is a 5 % damper on the foot (`FootDamping`).
- The pad-to-strut offsets (17–27 mm between strut lines at the foot) are lumped into one node.

## Not yet done

- The flight model's `Ground_Contact` still uses the sizer feet: 0.5 m radius, 17 MN/m, almost no damping. It should take this leg model (a MATLAB Function block) and the CAD feet (0.858 m radius, 27° off the body axes).
- Monte Carlo touchdown states could feed `landingSim` for load statistics.
