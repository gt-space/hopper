# Purchase requirements generated from Hopper FullFlow

Source config (relative to this script): ../../hopper_config.json
Run settings: dt=0.01 s; vent target=60 s; drain target=60 s; OFD fill is a CdA-versus-time trade study with no target fill time.

## Common pressure requirement
- Minimum working-pressure rating: 625 psia (610 psig), using 1.25x the configured 500 psia operating pressure.
- Confirm proof/burst ratings, pressure at the actual cryogenic temperature, port/thread rating, and system relief setpoint separately.

## OVNT (LOX ullage nitrogen vent)
- Effective flow area: at least 0.60 mm2 for the assumed 60 s vent target; configured candidate is 9.40 mm2 (water-equivalent Cv about 0.55).
- Fluid/service: dry nitrogen ullage in LOX service; materials and seals must be oxygen-clean and compatible with LOX-adjacent hardware. Minimum modeled temperature: 80 K (-193 C), including margin.
- Function: normally closed on-off vent valve, venting 500 psia initial tank pressure to ambient. Verify the required electrical actuation, default state, cycle life, and vent-line plumbing with the vehicle team.
- Initial-fill boiloff screen: at 5 psig tank control pressure, candidate OVNT capacity is compared with displaced ullage gas plus GOX boiloff. The configurable heat-load sweep spans 0 to 500 W; replace it with Hopper thermal data before using the result for selection.

## FVNT (fuel ullage nitrogen vent)
- Effective flow area: at least 0.49 mm2 for the assumed 60 s vent target; configured candidate is 15.00 mm2 (water-equivalent Cv about 0.88).
- Fluid/service: nitrogen ullage above RP-1 at 293 K. Verify wetted materials/seals for RP-1 vapor exposure and pressure rating.
- Function: normally closed on-off vent valve, venting 500 psia initial tank pressure to ambient. Confirm electrical actuation, default state, cycle life, and port geometry.

## OFD (LOX fill/drain)
- Effective flow area: use at least the worst-case 120.00 mm2 for the assumed 60 s drain target across the modeled pressure cases; configured candidate is 9.40 mm2 (water-equivalent Cv about 0.55).
- Fill trade study: fixed LFS inputs are 250 psig supply and 2.77 mm2 LFS path CdA. The OQD is a separate configurable 50.00 mm2 input. OFD CdA is swept and the resulting fill times are reported; no fill-time requirement or recommended valve size is imposed by this script.
- Fluid/service: liquid oxygen at 90 K. Require LOX compatibility, oxygen cleaning, cryogenic temperature rating, and a rated working pressure of at least 610 psig.
- Function: normally closed fill/drain valve. Confirm whether it also needs a separate ullage vent flow path; the FullFlow P&ID models OFD only on the liquid side.
- Conditional vent result: if OFD is plumbed to LOX ullage, its 9.40 mm2 CdA vents from 500 psia to 1 psig in 3.84 s. This is not a valid installed-path result until that ullage connection is confirmed.

## Model limits - do not use this as a qualification report
- No piping losses, valve travel, heat leak, two-phase flashing/cavitation, or back-pressure transient are modeled.
- The 8.5 in tank OD is used as an unverified effective ID for liquid-head calculation.
- The 60 s vent/drain targets and 48 min fill calibration are explicit sizing assumptions because the supplied Hopper configuration does not state them.
- Catalog Cv is only a water-flow comparison. Use vendor cryogenic/gas flow curves and tested valve CdA before selecting hardware.
