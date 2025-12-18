# -*- coding: utf-8 -*-
"""
Example usage of turbojet engine library

@author : Nathan_AZO

"""

import matplotlib.pyplot as plt
from lib import compute_engine_performance

# ----------------------
# Engine parameters
# ----------------------
dm = 210       # Air mass flow [kg/s]
beta = 5       # Bypass ratio
pi_fan = 1.5   # Fan compression ratio
pi_hpc_list = [n / 10 for n in range(10, 150, 10)]  # HPC compression ratio range
Tt4 = 1300     # Combustor exit temperature [K]
M0 = 0.85      # Flight Mach number
Ts0 = 280.65   # Static temperature at inlet [K]
Ps0 = 89876    # Static pressure at inlet [Pa]
M1 = 0.8
M2 = 0.5

# ----------------------
# Compute performance
# ----------------------
results = compute_engine_performance(dm, beta, pi_fan, pi_hpc_list, Tt4, M0, Ts0, Ps0, M1, M2)

# ----------------------
# Plot results
# ----------------------
x = pi_hpc_list
plt.figure(figsize=(12, 8))

# Subplot 1: Net and component thrusts
plt.subplot(2, 2, 1)
plt.plot(x, results["T_list"], label="Net Thrust")
plt.plot(x, results["T_fp_list"], label="Primary Flow Thrust")
plt.plot(x, results["T_fs_list"], label="Secondary Flow Thrust")
plt.title("Evolution of net and gross thrusts vs mass flow")
plt.legend()

# Subplot 2: Specific thrust and TSFC
plt.subplot(2, 2, 2)
plt.plot(x, results["ST_list"], label="Specific Thrust")
plt.plot(x, results["TSFC_list"], label="Specific Fuel Consumption")
plt.title("Evolution of specific thrust and TSFC vs mass flow")
plt.legend()

# Subplot 3: Powers
plt.subplot(2, 2, 3)
plt.plot(x, results["Wprop_list"], label="Propulsive Power")
plt.plot(x, results["Wcin_list"], label="Kinetic Power")
plt.plot(x, results["Wth_list"], label="Thermal Power")
plt.title("Evolution of powers vs mass flow")
plt.legend()

# Subplot 4: Efficiencies
plt.subplot(2, 2, 4)
plt.plot(x, results["nprop_list"], label="Propulsive Efficiency")
plt.plot(x, results["nth_list"], label="Thermal Efficiency")
plt.plot(x, results["nglob_list"], label="Global Efficiency")
plt.title("Evolution of efficiencies vs mass flow")
plt.legend()

plt.tight_layout()
plt.show()
