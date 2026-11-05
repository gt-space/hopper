"""Validate a complete saved feature-check run (nonzero regen heating required).

Usage: python tests/check_integration.py validation_outputs/burn_optimized.h5
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
from Hopper_Export import LoadTracks, _track, ExportActualPerformance
from Hopper_Report import summarize
from Hopper_Prop_System import heated_property_function, liquid_properties


def check(filename):
    time, tracks = LoadTracks(filename)
    get = lambda name: _track(tracks, name)
    dt = np.diff(time)
    assert np.all(dt > 0)
    pc = get("Chamber Pressure [psia]") * 6894.757293168
    rise = get("Regen Jacket Outlet Temperature [K]") - get("Regen Jacket Inlet Temperature [K]")
    assert np.max(rise) > 1.0, "Validation run must exercise nonzero heating."
    np.testing.assert_allclose(rise, np.interp(pc, [0, 2e6], [0, 60]), atol=1e-7)
    np.testing.assert_allclose(get("FINJ Inlet Temperature [K]"), get("Regen Jacket Outlet Temperature [K]"))
    assert np.max(get("Regen Jacket Heat Pickup [W]")) > 0
    for label in ("LOX", "Fuel"):
        mass = get(f"{label} Liquid Mass [kg]")
        flow = get(f"{label} Tank Outflow [kg/s]")
        np.testing.assert_allclose(mass[0] - mass[-1], np.sum(flow[1:] * dt), atol=2e-4)
    for name, opening_time, closing_time in [("OBB", .05, .05), ("FBB", .05, .05),
                                             ("OTV", .2, .2), ("FTV", .2, .2)]:
        opening = get(f"{name} Opening")
        assert np.min(opening) >= -1e-12 and np.max(opening) <= 1 + 1e-12
        change = np.diff(opening)
        assert np.all(change <= dt / opening_time + 1e-8)
        assert np.all(change >= -dt / closing_time - 1e-8)
    for name in ("COPV 1", "COPV 2", "COPV 3"):
        mass, flow = get(f"{name} Gas Mass [kg]"), get(f"{name} Outflow [kg/s]")
        # At discrete valve events, the saved connector rate is the endpoint
        # value, not the interval-average transfer. Check source direction and
        # conserved inventory over the connected gas network below.
        assert np.min(flow) >= -1e-12
        assert np.min(np.diff(mass)) < 0
        assert np.max(np.diff(mass)) < 5e-6
    np.testing.assert_allclose(get("COPV Gas Mass [kg]"),
                              sum(get(f"COPV {i} Gas Mass [kg]") for i in (1, 2, 3)), atol=1e-10)
    gas_mass = (get("COPV Gas Mass [kg]") + get("COPV Header Gas Mass [kg]")
                + get("LOX Ullage Gas Mass [kg]") + get("Fuel Ullage Gas Mass [kg]"))
    np.testing.assert_allclose(gas_mass, gas_mass[0], atol=3e-4)
    report = summarize(filename)
    assert report["copv_temperatures"][0]["violated"]
    assert report["copv_temperatures"][0]["limit_K"] == 294
    assert report["demonstrated_maximum_N"] is not None
    assert not report["invalid_inventory"]
    assert report["profile_completed"]
    # Check precomputed property accuracy at points between grid nodes.
    interpolation = heated_property_function("RP-1", 500*6894.757293168, 293, 353)
    for temperature in (293.3, 307.7, 332.9, 352.7):
        direct = liquid_properties("RP-1", 500*6894.757293168, temperature)
        approximate = interpolation(temperature)
        for name in ("density", "dynamic_viscosity", "enthalpy"):
            np.testing.assert_allclose(approximate[name], direct[name], rtol=2e-3)
    ExportActualPerformance(filename, str(Path(filename).parent / "validated_performance.csv"))
    print("Integration checks passed: regen, valve rates, inventory conservation, temperature limits, report, property accuracy, CSV export.")


if __name__ == "__main__":
    check(sys.argv[1])
