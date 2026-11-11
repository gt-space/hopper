# Purchase requirements generated from Hopper FullFlow

Source config (relative to this script): ../../hopper_config.json
Run settings: dt=0.01 s; vent target=60 s; drain target=60 s.

## Common pressure requirement
- Minimum working-pressure rating: 625 psia (610 psig), using 1.25x the configured 500 psia operating pressure.
- Confirm proof/burst ratings, pressure at the actual cryogenic temperature, port/thread rating, and system relief setpoint separately.

## OVNT (LOX ullage nitrogen vent)
- Effective flow area: at least 0.60 mm2 for the assumed 60 s vent target; configured candidate is 9.40 mm2 (water-equivalent Cv about 0.55).
- Fluid/service: dry nitrogen ullage in LOX service; materials and seals must be oxygen-clean and compatible with LOX-adjacent hardware. Minimum modeled temperature: 80 K (-193 C), including margin.
- Function: normally closed on-off vent valve, venting 500 psia initial tank pressure to ambient. Verify the required electrical actuation, default state, cycle life, and vent-line plumbing with the vehicle team.

## FVNT (fuel ullage nitrogen vent)
- Effective flow area: at least 0.49 mm2 for the assumed 60 s vent target; configured candidate is 15.00 mm2 (water-equivalent Cv about 0.88).
- Fluid/service: nitrogen ullage above RP-1 at 293 K. Verify wetted materials/seals for RP-1 vapor exposure and pressure rating.
- Function: normally closed on-off vent valve, venting 500 psia initial tank pressure to ambient. Confirm electrical actuation, default state, cycle life, and port geometry.

## OFD (LOX fill/drain)
- Effective flow area: use at least the worst-case 120.00 mm2 for the assumed 60 s drain target across the modeled pressure cases; configured candidate is 9.40 mm2 (water-equivalent Cv about 0.55).
- Fluid/service: liquid oxygen at 90 K. Require LOX compatibility, oxygen cleaning, cryogenic temperature rating, and a rated working pressure of at least 610 psig.
- Function: normally closed fill/drain valve. Confirm whether it also needs a separate ullage vent flow path; the FullFlow P&ID models OFD only on the liquid side.

## Model limits - do not use this as a qualification report
- No piping losses, valve travel, heat leak, two-phase flashing/cavitation, or back-pressure transient are modeled.
- The 8.5 in tank OD is used as an unverified effective ID for liquid-head calculation.
- The 60 s vent/drain targets are explicit sizing assumptions because the supplied configuration does not state these requirements.
- Catalog Cv is only a water-flow comparison. Use vendor cryogenic/gas flow curves and tested valve CdA before selecting hardware.
