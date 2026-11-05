"""Generate a box-and-arrow hardware diagram directly from the configured P&ID."""

from pathlib import Path
import textwrap


def hardware_diagram(output_dir="Hopper_Results"):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle, FancyArrowPatch
    import HOPPER_DATA as d

    path_length = max(len(d.LOX_PRESSURANT_LEG) + 1 + len(d.LOX_FEED),
                      len(d.FUEL_PRESSURANT_LEG) + 1 + len(d.FUEL_FEED))
    chamber_x = max(15.4, 1.45 * path_length + 0.9)
    fig, ax = plt.subplots(figsize=(max(26, chamber_x * 1.55), 10))
    ax.set_xlim(-0.9, max(chamber_x+1.4, len(d.COPVS)*1.65))
    ax.set_ylim(-2.3, 4.8)
    ax.axis("off")
    width, height = 1.12, 0.98

    def box(x, y, title, details):
        ax.add_patch(Rectangle((x-width/2, y-height/2), width, height,
                               facecolor="white", edgecolor="0.3", linewidth=1))
        title = "\n".join(textwrap.wrap(title, 19))
        ax.text(x, y + 0.30, title, ha="center", va="center", fontsize=8, weight="bold")
        ax.text(x, y - 0.11, "\n".join(details), ha="center", va="center", fontsize=7.5)

    def arrow(a, b):
        ax.add_patch(FancyArrowPatch(a, b, arrowstyle="->", mutation_scale=12,
                                     linewidth=1, color="0.35", connectionstyle="arc3"))

    def details(element):
        rows = [type(element).__name__]
        if hasattr(element, "inner_diameter"):
            rows.append(f"ID {element.inner_diameter*1000:.3g} mm")
            if getattr(element, "outer_diameter", None) is not None:
                rows.append(f"OD/wall {element.outer_diameter*1000:.3g}/{element.wall_thickness*1000:.3g} mm")
            rows.append(f"Length {element.length:g} m")
        elif hasattr(element, "cda_value"):
            rows.append(f"CdA {element.cda_value*1e6:.3g} mm2")
        elif hasattr(element, "max_cda"):
            rows.append(f"Max CdA {element.max_cda*1e6:.3g} mm2")
        if hasattr(element, "open_time"):
            rows.append(f"Open/close {element.open_time:g}/{element.close_time:g} s")
        return rows

    # Every bottle is a separate source, feeding one common header.
    count = len(d.COPVS)
    for i, bottle in enumerate(d.COPVS):
        x = i * 1.65
        limit = "off" if bottle.min_temperature is None else f"{bottle.min_temperature:g} K"
        box(x, 3.5, bottle.name, [d.PRESSURANT, f"Volume {bottle.volume*1000:g} L",
                                 f"Initial {bottle.pressure/6894.757:g} psia", f"Minimum T {limit}"])
        arrow((x, 3.5-height/2), (4.8-width/2, 2.55))
    box(4.8, 2.55, "Common header", ["Adiabatic gas storage", f"Volume {d.COPV_HEADER_VOLUME*1000:g} L"])

    for y, label, tank, leg, feed in [(1.1, "LOX", d.LOX_TANK, d.LOX_PRESSURANT_LEG, d.LOX_FEED),
                                     (-1.0, "Fuel", d.FUEL_TANK, d.FUEL_PRESSURANT_LEG, d.FUEL_FEED)]:
        xs = [i * 1.45 for i in range(len(leg) + 1 + len(feed))]
        tank_index = len(leg)
        # Route the header branch outside the component boxes.
        ax.plot([4.8-width/2, -0.78, -0.78], [2.55, 2.55, y], color="0.35", linewidth=1)
        arrow((-0.78, y), (xs[0]-width/2, y))
        entries = list(leg) + [tank] + list(feed)
        for i, (x, element) in enumerate(zip(xs, entries)):
            if i == tank_index:
                rows = [tank.propellant, f"Volume {tank.volume*1000:g} L", f"Initial mass {tank.liquid_mass:g} kg",
                        f"Set {tank.set_pressure/6894.757:g} psia"]
            else:
                rows = details(element)
            box(x, y, element.name, rows)
            if i:
                arrow((xs[i-1]+width/2, y), (x-width/2, y))
        arrow((xs[-1]+width/2, y), (chamber_x-width/2, 0.05))
    box(chamber_x, 0.05, "Main chamber / nozzle", [f"O/F target {d.MIXTURE_RATIO:g}",
                f"Throat {d.THROAT_AREA*1e6:.3g} mm2", f"Area ratio {d.EXPANSION_RATIO:.3g}"])
    ax.text(-0.55, 4.55, "Hopper propulsion system", fontsize=16, weight="bold")
    ax.text(-0.55, 4.30, "Sources and components from HOPPER_DATA.py; dimensions and timings are configured inputs.", fontsize=10)
    ax.text(-0.55, -2.05, "Arrows show connectivity. Layout is schematic; component position does not represent physical distance.", fontsize=9)
    fig.tight_layout()
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    for extension in ("png", "svg"):
        fig.savefig(output / f"component_diagram.{extension}", dpi=150)
    plt.close(fig)
    return output / "component_diagram.svg"


if __name__ == "__main__":
    import HOPPER_DATA as data
    print(hardware_diagram(data.PROP_RESULTS_DIRECTORY))
