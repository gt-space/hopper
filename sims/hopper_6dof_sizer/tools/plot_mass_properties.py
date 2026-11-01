"""Hopper mass properties over the nominal flight: CAD dry block vs the 116.86 kg stand-in.

Python port of the sim's mass model so it runs without MATLAB:
  - vehicle_mass_props.m: dry block + LOx/RP-1 liquid columns on the tank floors
  - cg_moi_init.m (empty tanks, inertia scaled to the dry mass): the pre-CAD
    dry block, as mass_properties_setup.m built it
  - Trajectory2.m: thrust reference m(g + z''), clipped and low-passed

Propellant flow is proportional to the thrust reference. Its scale is fixed
by the nominal baseline (tests/baselines/nominal_v2.mat: 28.39 s flight,
3.520 kg LOx and 3.016 kg RP-1 left at touchdown) for the stand-in vehicle;
the CAD vehicle uses the same Isp, so its flow follows its own mass.

Values mirror vehicles/hopper/hopper_config.m (CAD 2026-10-29) and
inputs/mission_inputs.m; update them here if those change.

    python tools/plot_mass_properties.py [out_dir]   (default outputs/figures)
"""
import sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

G = 9.81

# ---- Vehicle (hopper_config.m) ----
TANK_R = 4.25 * 0.0254
OX = dict(h=29.18 * 0.0254, cg=1.854, wall=0.003175, mass=4.6940418420227319, rho=1141.0, load=16.0)
FU = dict(h=18.91 * 0.0254, cg=1.0668, wall=0.003175, mass=5.6082853202836009, rho=800.0, load=13.0)
for t in (OX, FU):
    t["floor"] = t["cg"] - t["h"] / 2

# CAD dry block (body frame: x up the thrust axis from the footpad bottoms)
TO_BODY = np.array([[0, 1, 0], [0, 0, 1], [1, 0, 0]], float)
CAD_MASS = 114.894
CAD_CG = np.array([-0.003159, -0.169015, 0.003797])
CAD_I = np.array([[51.5454, 0.2508, -0.0406], [0.2508, 12.6948, 0.1037], [-0.0406, 0.1037, 52.3340]])
Y_BASE = -1.208957
cg_b = TO_BODY @ (CAD_CG - np.array([0, Y_BASE, 0]))
CAD = dict(mass=CAD_MASS, cg=np.array([cg_b[0], 0, 0]), I=np.diag(np.diag(TO_BODY @ CAD_I @ TO_BODY.T)),
           lat=cg_b[1:])

# ---- Nominal trajectory (mission_inputs.m) and baseline ----
APEX, T_UP, T_HOVER, T_DOWN, T_TERM, Z_TD, V_TD = 52, 12, 2, 13, 2.5, 0.4, 0.3
BURN_MASS_REF = 21.0
THRUST_NOM, THROTTLE = 500 * 4.44822, (0.4, 1.1)
T_FLIGHT = 28.39
OX_END, FU_END = 3.519973779765257, 3.015979023814346
STANDIN_DRY = 116.86


def standin_dry_block(dry_mass):
    """cg_moi_init with empty tanks, inertia scaled to dry_mass (mass_properties_setup.m)."""
    al = 2700.0
    eng = (37.83, 0.3556)
    others = [(34.26, 1.193), (11.52, 1.4224), (28.41, 1.0668), (15.0, 2.0)]  # structures, avionics, fluids, payload
    copv_m, copv_r, copv_t, copv_h, copv_x, copv_y = 11 + 6, 0.18, 13.21e-3, 0.56, 1.21, 0.23
    it_h = TANK_R
    caps = []
    for t in (OX, FU):
        mc = al * np.pi * TANK_R**2 * t["wall"]
        caps += [(mc, t["cg"] - t["h"] / 2 - t["wall"] / 2, t["wall"]), (mc, t["cg"] + t["h"] / 2 + t["wall"] / 2, t["wall"])]
    pts = [eng, (FU["mass"], FU["cg"]), (OX["mass"], OX["cg"]), (copv_m, copv_x)] + others + [c[:2] for c in caps]
    m_comp = sum(m for m, _ in pts)
    cg = sum(m * x for m, x in pts) / m_comp

    s_izz = lambda m, r: 0.5 * m * r**2
    s_ixx = lambda m, r, h: m * (3 * r**2 + h**2) / 12
    h_izz = lambda m, ro, ri: 0.5 * m * (ro**2 + ri**2)
    h_ixx = lambda m, ro, ri, h: m * (3 * (ro**2 + ri**2) + h**2) / 12
    eng_r, eng_h = 4 * 0.0254, 0.3048
    izz = s_izz(eng[0], eng_r)
    ixx = s_ixx(eng[0], eng_r, eng_h) + eng[0] * (eng[1] - cg) ** 2
    for t in (FU, OX):
        ri = TANK_R - t["wall"]
        izz += h_izz(t["mass"], TANK_R, ri)
        ixx += h_ixx(t["mass"], TANK_R, ri, t["h"]) + t["mass"] * (t["cg"] - cg) ** 2
    ixx_copv = h_ixx(copv_m, copv_r, copv_r - copv_t, copv_h) + copv_m * (copv_x - cg) ** 2
    izz += h_izz(copv_m, copv_r, copv_r - copv_t)
    ixx += ixx_copv
    for m, x in others:
        izz += s_izz(m, TANK_R)
        ixx += s_ixx(m, TANK_R, it_h) + m * (x - cg) ** 2
    for m, x, w in caps:
        izz += s_izz(m, TANK_R)
        ixx += s_ixx(m, TANK_R, w) + m * (x - cg) ** 2
    iyy = ixx + 2 * ixx_copv * copv_y**2   # as cg_moi_init.m has it
    I = np.diag([izz, iyy, ixx]) * (dry_mass / m_comp)
    return dict(mass=dry_mass, cg=np.array([cg, 0, 0]), I=I, lat=np.zeros(2))


def vehicle_mass_props(dry, m_ox, m_fu):
    """vehicle_mass_props.m: CG [3] and inertia about it [3x3], body axes."""
    A = np.pi * TANK_R**2
    h_ox, h_fu = m_ox / (OX["rho"] * A), m_fu / (FU["rho"] * A)
    r_ox = np.array([OX["floor"] + h_ox / 2, 0, 0])
    r_fu = np.array([FU["floor"] + h_fu / 2, 0, 0])
    M = dry["mass"] + m_ox + m_fu
    cg = (dry["mass"] * dry["cg"] + m_ox * r_ox + m_fu * r_fu) / M
    pa = lambda m, d: m * ((d @ d) * np.eye(3) - np.outer(d, d))
    col = lambda m, h: np.diag([m * TANK_R**2 / 2, m * (3 * TANK_R**2 + h**2) / 12, m * (3 * TANK_R**2 + h**2) / 12])
    I = (dry["I"] + pa(dry["mass"], dry["cg"] - cg) + col(m_ox, h_ox) + pa(m_ox, r_ox - cg)
         + col(m_fu, h_fu) + pa(m_fu, r_fu - cg))
    return cg, I


def thrust_reference(m0, dt=0.01):
    """Trajectory2.m thrust reference for a vehicle of wet mass m0."""
    t_end = T_UP + T_HOVER + T_DOWN + T_TERM
    t = np.arange(0, t_end + dt / 2, dt)
    m = m0 - BURN_MASS_REF * t / t_end
    acc = np.zeros_like(t)
    up = t < T_UP
    s = t[up] / T_UP
    acc[up] = APEX * (60 * s - 180 * s**2 + 120 * s**3) / T_UP**2
    dn = (t >= T_UP + T_HOVER) & (t < T_UP + T_HOVER + T_DOWN)
    s = (t[dn] - T_UP - T_HOVER) / T_DOWN
    acc[dn] = ((Z_TD - APEX) * (60 * s - 180 * s**2 + 120 * s**3)
               - V_TD * T_DOWN * (-24 * s + 84 * s**2 - 60 * s**3)) / T_DOWN**2
    T = np.clip(m * (G + acc), THRUST_NOM * THROTTLE[0], THRUST_NOM * THROTTLE[1])
    a = dt / (0.2 + dt)
    for k in range(1, len(T)):
        T[k] = T[k - 1] + a * (T[k] - T[k - 1])
    return t, T


def history(dry, mdot_per_N):
    ox_frac = (OX["load"] - OX_END) / ((OX["load"] - OX_END) + (FU["load"] - FU_END))
    t, T = thrust_reference(dry["mass"] + OX["load"] + FU["load"])
    keep = t <= T_FLIGHT + 1e-9
    t, T = t[keep], T[keep]
    burned = np.concatenate([[0], np.cumsum((T[1:] + T[:-1]) / 2 * np.diff(t))]) * mdot_per_N
    m_ox = OX["load"] - ox_frac * burned
    m_fu = FU["load"] - (1 - ox_frac) * burned
    out = dict(t=t, mass=dry["mass"] + m_ox + m_fu, cg=[], lat=[], I=[])
    for mo, mf in zip(m_ox, m_fu):
        cg, I = vehicle_mass_props(dry, mo, mf)
        out["cg"].append(cg[0])
        out["I"].append(np.diag(I))
        out["lat"].append(dry["lat"] * dry["mass"] / (dry["mass"] + mo + mf))
    for k in ("cg", "lat", "I"):
        out[k] = np.array(out[k])
    return out


def main(out_dir):
    out_dir.mkdir(parents=True, exist_ok=True)
    old = standin_dry_block(STANDIN_DRY)
    # flow scale from the baseline (stand-in vehicle): burned propellant at touchdown
    t, T = thrust_reference(STANDIN_DRY + OX["load"] + FU["load"])
    keep = t <= T_FLIGHT + 1e-9
    impulse = np.trapezoid(T[keep], t[keep])
    mdot_per_N = ((OX["load"] - OX_END) + (FU["load"] - FU_END)) / impulse   # 1 / (Isp g0)
    new, ref = history(CAD, mdot_per_N), history(old, mdot_per_N)

    BLUE, ORANGE = "#2a78d6", "#eb6834"
    INK, INK2, GRID, SURF = "#0b0b0b", "#52514e", "#e6e5e0", "#fcfcfb"
    plt.rcParams.update({"font.size": 10, "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2,
                         "ytick.color": INK2, "axes.spines.top": False, "axes.spines.right": False,
                         "figure.facecolor": SURF, "axes.facecolor": SURF, "savefig.facecolor": SURF})

    # One figure, CAD dry block only: mass, CG height, pitch & yaw, roll
    f, axes = plt.subplots(4, 1, figsize=(8, 10.5), sharex=False)
    f.suptitle("CAD mass properties, nominal flight", fontsize=13, fontweight="bold", color=INK)
    t = new["t"]

    def panel(ax, title, ylabel):
        ax.set_title(title, fontsize=11, color=INK)
        ax.set_ylabel(ylabel)
        ax.grid(True, color=GRID, lw=0.8)
        ax.set_axisbelow(True)
        ax.set_xlim(0, 31)
        ax.margins(y=0.12)

    def end_label(ax, y, text, dy=0):
        ax.annotate(text, (t[-1], y[-1]), xytext=(6, dy), textcoords="offset points",
                    va="center", fontsize=9, color=INK2)

    ax = axes[0]
    panel(ax, "Total mass", "Mass [kg]")
    ax.plot(t, new["mass"], color=BLUE, lw=2)
    end_label(ax, new["mass"], f"{new['mass'][-1]:.1f} kg")

    ax = axes[1]
    panel(ax, "CG height", "CG above base [m]")
    ax.plot(t, new["cg"], color=BLUE, lw=2)
    end_label(ax, new["cg"], f"{new['cg'][-1]:.3f} m")

    ax = axes[2]
    panel(ax, "Pitch and yaw inertia", "Inertia [kg·m²]")
    ax.plot(t, new["I"][:, 1], color=BLUE, lw=2, label="Pitch")
    ax.plot(t, new["I"][:, 2], color=ORANGE, lw=2, ls="--", label="Yaw")
    end_label(ax, new["I"][:, 1], f"{new['I'][-1, 1]:.2f}", 6)
    end_label(ax, new["I"][:, 2], f"{new['I'][-1, 2]:.2f}", -6)
    ax.legend(frameon=False, loc="upper right", labelcolor=INK)

    ax = axes[3]
    panel(ax, "Roll inertia", "Inertia [kg·m²]")
    ax.plot(t, new["I"][:, 0], color=BLUE, lw=2)
    end_label(ax, new["I"][:, 0], f"{new['I'][-1, 0]:.3f}")
    ax.set_xlabel("Time [s]")

    f.tight_layout()
    f.savefig(out_dir / "cad_mass_props_vs_time.png", dpi=150)
    plt.close(f)

    for name, h in (("CAD", new), ("stand-in", ref)):
        print(f"{name:9s} mass {h['mass'][0]:.3f}->{h['mass'][-1]:.3f} kg  cg {h['cg'][0]:.4f}->{h['cg'][-1]:.4f} m  "
              f"I0 [{', '.join(f'{v:.3f}' for v in h['I'][0])}]  Iend [{', '.join(f'{v:.3f}' for v in h['I'][-1])}]")
    print("lateral CG at t=0 [mm]:", np.round(new["lat"][0] * 1e3, 2), " Isp-equivalent:", round(1 / (mdot_per_N * 9.80665), 1), "s")


if __name__ == "__main__":
    main(Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).resolve().parents[1] / "outputs" / "figures")
