from rocketcea.cea_obj import CEA_Obj
import numpy as np

# -------------------------
# USER INPUTS
# -------------------------
ox = 'LOX'
fuel = 'RP1'

MR = 1.25         # constant mixture ratio
Pc_min = 100      # psia
Pc_max = 350
Pc_step = 10

# -------------------------
# CREATE CEA OBJECT
# -------------------------
cea = CEA_Obj(oxName=ox, fuelName=fuel)

# Pressure vector
Pc_vec = np.arange(Pc_min, Pc_max + Pc_step, Pc_step)

# Compute cstar vector (convert ft/s -> m/s)
Cstar_vec = []

for Pc in Pc_vec:
    cstar = cea.get_Cstar(
        Pc=float(Pc),
        MR=MR
    ) * 0.3048

    Cstar_vec.append(cstar)

# -------------------------
# PRINT MATLAB STYLE VECTORS
# -------------------------
Pc_string = "[" + " ".join(f"{x:.6g}" for x in Pc_vec) + "]"
Cstar_string = "[" + " ".join(f"{x:.6f}" for x in Cstar_vec) + "]"

print("Pc =")
print(Pc_string)

print("\nCstar_mps =")
print(Cstar_string)