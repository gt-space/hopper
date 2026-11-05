"""Run a short, isolated coupled feature check without editing user inputs.

python tests/run_feature_validation.py
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))


def main():
    import matplotlib
    matplotlib.use("Agg")
    import HOPPER_DATA as data
    import Hopper_Prop_System as simulation
    from Hopper_Elements import Regen
    from Hopper_Thrust_Targeter import HopperThrustTargets
    from Hopper_Report import export_summary
    from Hopper_Plots import HopperRun, make_plots
    from check_integration import check

    output = Path("validation_outputs")
    output.mkdir(exist_ok=True)
    data.PROP_THRUST_CURVE_NAME = "feature-validation"
    for element in simulation.FUEL_FEED:
        if isinstance(element, Regen):
            element.delta_t_curve = ((0., 0.), (2e6, 60.))
    simulation.COPVS[0].min_temperature = 294.  # deliberately above initial temperature
    targets = HopperThrustTargets({0.: 500., .5: 700.})
    filename = str(output / "burn_optimized.h5")
    print("Running 1.5 s startup/burn/shutdown feature check at <=5 ms integration steps...")
    simulation.HopperPropSystem(targets, filename=filename, save_dt=.005, statistics=False)
    export_summary(filename, output / "burn")
    check(filename)
    paths = make_plots(HopperRun(filename), "all", output / "plots", csv_out=True)
    print(f"Rendered {len(paths)//2} plots with CSVs in {output / 'plots'}")


if __name__ == "__main__":
    main()
