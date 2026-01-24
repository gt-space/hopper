import CoolProp.CoolProp as CP
import numpy as np

# --- CONFIGURATION ---
burn_duration = 25.0       # s
mdot_lox = 0.777           # kg/s
mdot_fuel = 0.388          # kg/s

# Pressures
P_tank_target = 500.0      # psi (Regulated pressure)
P_eol_margin = 0.00        # 20% margin on tank pressure
P_copv_start = 4500.0      # psi
P_regulator_drop = 50.0    # psi (Min pressure drop across regulator)

# Constraints
T_min_copv = 233.0         # K (Minimum allowable EOL Temp for COPV material)

# Geometry / Physics
ullage_frac = 0.15         # 15% Ullage
collapse_factor_lox = 2.0  # LOX Tank Factor (applied to COLD gas density)
collapse_factor_fuel = 1.0 # Fuel Tank Factor

# Temperatures
T_lox = 90.0               # K
T_fuel = 293.0             # K
T_copv_init = 293.0        # K

# Units
PSI_TO_PA = 6894.75729

# --- CALCULATIONS ---

P_tank_eol_psi = P_tank_target * (1 + P_eol_margin)
P_copv_cutoff_psi = P_tank_eol_psi + P_regulator_drop

print(f"--- DESIGN CONSTRAINTS ---")
print(f"Prop Tank Target:   {P_tank_eol_psi:.0f} psi")
print(f"Regulator Cutoff:   {P_copv_cutoff_psi:.0f} psi")
print(f"Min Gas Supply T:   {T_min_copv:.0f} K (Used for density sizing)")
print("-" * 30)

def size_regulated_leg(fluid_name, mdot, T_fluid, collapse_factor):
    # 1. Propellant Volume
    m_prop = mdot * burn_duration
    rho_prop = CP.PropsSI('D', 'T', T_fluid, 'P', P_tank_eol_psi*PSI_TO_PA, fluid_name)
    V_prop = m_prop / rho_prop
    
    # 2. Tank Volume
    V_total = V_prop / (1 - ullage_frac)
    V_ullage = V_total * ullage_frac
    
    # 3. Gas Mass Required (UPDATED LOGIC)
    # Instead of assuming gas is 293K, we assume the ullage is filled with
    # gas at the MINIMUM supply temp (200K). This accounts for supply cooling.
    
    # Base Density (Cold Gas)
    rho_gn2_cold = CP.PropsSI('D', 'P', P_tank_eol_psi*PSI_TO_PA, 'T', T_min_copv, 'Nitrogen')
    
    # Apply Collapse Factor to the COLD density
    m_gas_final = V_total * rho_gn2_cold * collapse_factor
    
    # Initial mass (usually negligible or warm, but let's be conservative and subtract warm mass)
    rho_gn2_warm_init = CP.PropsSI('D', 'P', P_tank_eol_psi*PSI_TO_PA, 'T', T_copv_init, 'Nitrogen')
    m_gas_initial = V_ullage * rho_gn2_warm_init
    
    m_gas_transfer = max(0, m_gas_final - m_gas_initial)
    
    return {
        "V_total": V_total,
        "m_gas_transfer": m_gas_transfer,
        "m_prop": m_prop,
        "rho_base": rho_gn2_cold
    }

# --- SIZE PROPELLANT TANKS ---
lox = size_regulated_leg("Oxygen", mdot_lox, T_lox, collapse_factor_lox)
fuel = size_regulated_leg("n-Dodecane", mdot_fuel, T_fuel, collapse_factor_fuel)

# --- SIZE COPV (With Temp Constraint) ---
total_gas_mass_needed = lox['m_gas_transfer'] + fuel['m_gas_transfer']

# State 1: Start
rho_copv_1 = CP.PropsSI('D', 'P', P_copv_start*PSI_TO_PA, 'T', T_copv_init, 'Nitrogen')
s_copv_1   = CP.PropsSI('S', 'P', P_copv_start*PSI_TO_PA, 'T', T_copv_init, 'Nitrogen')

# State 2: Determine Residual Limit (Same Temp Logic as before)
# Check Isentropic Temp at Regulator Cutoff
T_at_cutoff = CP.PropsSI('T', 'P', P_copv_cutoff_psi*PSI_TO_PA, 'S', s_copv_1, 'Nitrogen')

if T_at_cutoff < T_min_copv:
    # Limit by Temperature
    limit_reason = "Temperature (200K)"
    T_residual = T_min_copv
    
    # Find P at T=200K isentropically
    p_low, p_high = P_copv_cutoff_psi*PSI_TO_PA, P_copv_start*PSI_TO_PA
    P_residual_pa = p_low
    
    for _ in range(20):
        p_guess = 0.5 * (p_low + p_high)
        s_guess = CP.PropsSI('S', 'P', p_guess, 'T', T_min_copv, 'Nitrogen')
        if abs(s_guess - s_copv_1) < 1.0:
            P_residual_pa = p_guess
            break
        elif s_guess > s_copv_1:
            p_low = p_guess
        else:
            p_high = p_guess
else:
    # Limit by Regulator
    limit_reason = "Regulator Pressure"
    P_residual_pa = P_copv_cutoff_psi * PSI_TO_PA
    T_residual = T_at_cutoff

# Residual Density
rho_residual = CP.PropsSI('D', 'P', P_residual_pa, 'T', T_residual, 'Nitrogen')

# Sizing
delta_rho = rho_copv_1 - rho_residual
V_copv_required = total_gas_mass_needed / delta_rho

# --- OUTPUT ---
def print_res(name, data):
    print(f"\n{name} TANK:")
    print(f"  Mass:         {data['m_prop']:.2f} kg")
    print(f"  Volume:       {data['V_total']*1000:.2f} L")
    print(f"  GN2 Needed:   {data['m_gas_transfer']:.3f} kg")
    print(f"  Base Density: {data['rho_base']:.2f} kg/m3 (at {T_min_copv}K)")

print_res("LOX", lox)
print_res("FUEL", fuel)

print(f"\nSUPPLY COPV:")
print(f"  Gas Mass Req:   {total_gas_mass_needed:.2f} kg")
print(f"  Volume Required:  {V_copv_required*1000:.2f} Liters")
print(f"  Residual Press:   {P_residual_pa/PSI_TO_PA:.0f} psi")
print(f"  Residual Temp:    {T_residual:.1f} K")
print(f"  Total Gas Mass:   {CP.PropsSI('D', 'P', 4500*PSI_TO_PA, 'T', 293, 'Nitrogen')*V_copv_required:.2f} kg")
