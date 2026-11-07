"""Compact hardware schematic generated from the loaded JSON configuration."""
from pathlib import Path
import math
import sys
import textwrap

if not __package__:   # run as a script (python Backend/Hopper_Diagram.py): make HOPPER_DATA and the Backend
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))   # package importable


def hardware_diagram(output_dir="Outputs"):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Polygon, Ellipse
    import HOPPER_DATA as d

    branches = [("OXIDIZER", "#177f9c", d.LOX_PRESSURANT_LEG + [d.LOX_TANK] + d.LOX_FEED),
                ("FUEL", "#ba622a", d.FUEL_PRESSURANT_LEG + [d.FUEL_TANK] + d.FUEL_FEED)]
    rows = max(math.ceil(len(items)/3) for _, _, items in branches)
    top = rows*1.95 + 3.8 + math.ceil(len(d.COPVS)/6)*1.6
    fig, ax = plt.subplots(figsize=(17, top*.87), facecolor="#f5f7fb")
    ax.set(xlim=(-.7, 17.7), ylim=(-1.05, top+.6))
    ax.set_aspect("equal")
    ax.axis("off")
    width, height = 2.35, 1.48
    ink = "#233449"

    def arrow(points, color=ink):
        if len(points) > 2:
            ax.plot(*zip(*points[:-1]), color=color, lw=1.5, solid_capstyle="round")
        ax.add_patch(FancyArrowPatch(points[-2], points[-1], arrowstyle="-|>",
                                    mutation_scale=11, linewidth=1.5, color=color))

    def card(x, y, title, details, kind, color):
        ax.add_patch(FancyBboxPatch((x-width/2, y-height/2), width, height,
            boxstyle="round,pad=0.04,rounding_size=0.12", facecolor="white", edgecolor="#ced7e2", lw=.9))
        sx, sy = x-width/2+.23, y+.48
        if kind in ("Tank", "COPVSource"):
            ax.add_patch(FancyBboxPatch((sx-.10, sy-.14), .2, .28,
                boxstyle="round,pad=0,rounding_size=0.07", ec=color, fc="#e6eef7", lw=1.3))
        elif "Valve" in kind:
            ax.add_patch(Polygon([(sx-.13,sy-.10),(sx+.13,sy+.10),(sx+.13,sy-.10),(sx-.13,sy+.10)], ec=color, fc=color, lw=1))
        elif "Line" in kind:
            for dy in (-.045, .045):
                ax.plot([sx-.13,sx+.13],[sy+dy,sy+dy],color=color,lw=1.6)
        else:
            ax.add_patch(Ellipse((sx,sy),.23,.23,ec=color,fc="white",lw=1.4))
            ax.plot([sx,sx],[sy-.13,sy+.13],color=color,lw=1.2)
        ax.text(x-width/2+.45, y+.49, "\n".join(textwrap.wrap(title,23)),
                color=ink, fontsize=8.2, weight="bold", va="center")
        ax.text(x-width/2+.16, y+.12, "\n".join(details), color="#43556b",
                fontsize=7.6, va="top", linespacing=1.4)

    def details(e):
        if hasattr(e,"liquid_mass"):
            return [f"{e.propellant} | {e.volume*1000:g} L", f"Initial mass {e.liquid_mass:g} kg",
                    f"Setpoint {e.set_pressure/6894.757293168:g} psia"]
        if hasattr(e,"inner_diameter"):
            od, wall = getattr(e,"outer_diameter",None), getattr(e,"wall_thickness",None)
            return [f'OD {od/.0254:.4g} in | wall {wall/.0254:.4g} in' if od is not None else 'OD / wall: unspecified',
                    f'ID {e.inner_diameter/.0254:.4g} in', f'Length {e.length/.0254:.4g} in']
        if hasattr(e,"max_cda"):
            lines = [f"Max CdA {e.max_cda*1e6:.4g} mm²"]
        elif hasattr(e,"cda_value"):
            lines = [f"CdA {e.cda_value*1e6:.4g} mm²"]
        else:
            lines = [f"Throat {e.throat_area/.0254**2:.4g} in²"]
        if hasattr(e,"open_time"):
            lines += [f"Open {e.open_time:g} s / close {e.close_time:g} s"]
        if getattr(e,"length",None):
            lines += [f"Passage {e.length/.0254:.4g} in"]
        if hasattr(e,"delta_t_curve"):
            lines += ["Heating: configured ΔT(Pc)"]
        return lines

    ax.text(-.35,top+.22,"HOPPER  /  PROPULSION NETWORK",fontsize=17,weight="bold",color=ink)
    ax.text(-.35,top-.13,f"{d.CONFIG_PATH.name}  •  {d.PROP_THRUST_CURVE_NAME} demand  •  tube dimensions in inches",
            fontsize=9,color="#607287")
    header_y = rows*1.95+2.1
    for i,b in enumerate(d.COPVS):
        col,row = i%6,i//6
        x,y = 1.15+col*2.85,top-1.1-row*1.6
        limit = "disabled" if b.min_temperature is None else f"{b.min_temperature:g} K"
        card(x,y,b.name,[f"{d.PRESSURANT} | {b.volume*1000:g} L",
             f"{b.pressure/6894.757293168:g} psia / {b.temperature:g} K",f"Min T {limit}"],"COPVSource","#546b88")
        rail=x+width/2+.14
        arrow([(x+width/2,y),(rail,y),(rail,header_y+1),(8.5,header_y+1),(8.5,header_y+height/2)],"#75869b")
    card(8.5,header_y,"Common header",[f"Volume {d.COPV_HEADER_VOLUME*1000:g} L","Adiabatic gas storage"],"Tank","#546b88")
    first_y = rows*1.95+.25
    chamber_y = -.05
    for side,(label,color,items) in enumerate(branches):
        base=.95+side*9.05
        ax.text(base-1.05,first_y+1.04,label,fontsize=10,weight="bold",color=color)
        positions=[]
        for i,e in enumerate(items):
            row,col=divmod(i,3)
            col=col if row%2==0 else 2-col
            x,y=base+col*2.85,first_y-row*1.95
            positions.append((x,y))
            card(x,y,e.name,details(e),type(e).__name__,color)
            if i:
                px,py=positions[-2]
                if py==y:
                    sign=1 if x>px else -1
                    arrow([(px+sign*width/2,py),(x-sign*width/2,y)],color)
                else:
                    arrow([(px,py-height/2),(x,y+height/2)],color)
        x,y=positions[0]
        arrow([(8.5,header_y-height/2),(8.5,first_y+1.0),(x,first_y+1.0),(x,y+height/2)],color)
        x,y=positions[-1]
        sign=-1 if side==0 else 1
        arrow([(x,y-height/2),(x,chamber_y),(8.5+sign*width/2,chamber_y)],color)
    card(8.5,chamber_y,"Chamber / nozzle",[f"Throat area {d.THROAT_AREA/.0254**2:.5g} in²",
        f"Throat diameter {math.sqrt(4*d.THROAT_AREA/math.pi)/.0254:.5g} in",
        f"O/F {d.MIXTURE_RATIO:g} | area ratio {d.EXPANSION_RATIO:.4g}"],"Nozzle",ink)
    fig.tight_layout(pad=.8)
    output=Path(output_dir)
    output.mkdir(parents=True,exist_ok=True)
    for extension in ("png","svg"):
        fig.savefig(output/f"component_diagram.{extension}",dpi=150,facecolor=fig.get_facecolor())
    plt.close(fig)
    return output/"component_diagram.svg"


if __name__ == "__main__":
    import HOPPER_DATA as data
    print(hardware_diagram(data.PROP_RESULTS_DIRECTORY))
