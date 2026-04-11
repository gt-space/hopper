from rocketcea.cea_obj import CEA_Obj
import numpy as np

# -------------------------
# USER INPUTS
# -------------------------
ox = 'LOX'
fuel = 'RP1'

MR = 1.25
eps = 2.931          # nozzle expansion ratio Ae/At
Pamb = 14.7        # psia

Pc_min = 100       # psia
Pc_max = 350
Pc_step = 10

# -------------------------
# CREATE CEA OBJECT
# -------------------------
cea = CEA_Obj(oxName=ox, fuelName=fuel)

# Chamber pressure vector
Pc_vec = np.arange(Pc_min, Pc_max + Pc_step, Pc_step)

# Compute Cf vector
Cf_vec = []

for Pc in Pc_vec:
    CFcea, CFamb, mode = cea.get_PambCf(
        Pc=float(Pc),
        MR=MR,
        eps=eps,
        Pamb=Pamb
    )
    Cf_vec.append(CFamb)

# -------------------------
# PRINT MATLAB STYLE VECTORS
# -------------------------
Pc_string = "[" + " ".join(f"{x:.6g}" for x in Pc_vec) + "]"
Cf_string = "[" + " ".join(f"{x:.6f}" for x in Cf_vec) + "]"

print("Pc =")
print(Pc_string)

print("\nCf =")
print(Cf_string)