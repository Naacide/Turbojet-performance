# -*- coding: utf-8 -*-
"""
Turbojet engine performance calculator

This library computes performance parameters of a simple turbojet engine
with bypass (turbo-fan), using standard thermodynamic assumptions.

@author : Nathan_AZO
"""

import numpy as np

# ----------------------------------------------------------------------
# Engine constants
# ----------------------------------------------------------------------
gamma = 1.4  # Ratio of specific heats
R = 287      # Gas constant [J/kg/K]
PCI = 48e6   # Lower heating value of fuel [J/kg]
cp = 1004    # Specific heat at constant pressure [J/kg/K]


def compute_engine_performance(dm, beta, pi_fan, pi_hpc_list, Tt4, M0, Ts0, Ps0, M1, M2):
    """
    Computes the performance of a turbojet engine over a range of high-pressure compressor ratios.

    Parameters:
        dm (float): Air mass flow rate [kg/s]
        beta (float): Bypass ratio
        pi_fan (float): Fan compression ratio
        pi_hpc_list (list of float): List of high-pressure compressor ratios
        Tt4 (float): Combustor exit temperature [K]
        M0 (float): Flight Mach number
        Ts0 (float): Static temperature at inlet [K]
        Ps0 (float): Static pressure at inlet [Pa]
        M1, M2 (float): Mach numbers at station 1 and 2

    Returns:
        dict: Dictionary containing lists of performance parameters for each HPC ratio
    """
    # Initialize lists
    results = {
        "T_list": [],
        "T_fp_list": [],
        "T_fs_list": [],
        "ST_list": [],
        "TSFC_list": [],
        "Wprop_list": [],
        "Wcin_list": [],
        "Wth_list": [],
        "nprop_list": [],
        "nth_list": [],
        "nglob_list": [],
        "A1_list": [],
        "A2_list": [],
        "A8_list": [],
        "A18_list": []
    }

    for pi_hpc in pi_hpc_list:
        # ----------------------
        # Inlet conditions
        # ----------------------
        Tt2 = Ts0 * (1 + (gamma - 1) / 2 * M0 ** 2)
        Pt2 = Ps0 * (1 + (gamma - 1) / 2 * M0 ** 2) ** (gamma / (gamma - 1))

        # ----------------------
        # Fan exit
        # ----------------------
        Pt13 = Pt2 * pi_fan
        Tt13 = Tt2 * pi_fan ** ((gamma - 1) / gamma)
        dm_fp = dm / (1 + beta)   # primary flow
        dm_fs = dm - dm_fp        # secondary flow

        # ----------------------
        # High-pressure compressor exit
        # ----------------------
        Pt3 = Pt13 * pi_hpc
        Tt3 = Tt13 * pi_hpc ** ((gamma - 1) / gamma)

        # ----------------------
        # Combustor exit
        # ----------------------
        Pt4 = Pt3
        # Tt4 is given
        dm_f = (dm_fp * cp * (Tt4 - Tt3)) / (PCI - Tt4)
        f = dm_f / dm_fp  # fuel/air ratio

        # ----------------------
        # High-pressure turbine exit
        # ----------------------
        Tt42 = Tt4 - ((Tt3 - Tt13) / (1 + f))
        tau_hpt = Tt42 / Tt4
        Pt42 = Pt4 * tau_hpt ** (gamma / (gamma - 1))

        # ----------------------
        # Low-pressure turbine exit
        # ----------------------
        Tt5 = Tt42 - ((1 + beta) * (Tt13 - Tt2)) / (1 + f)
        tau_lpt = Tt5 / Tt42
        Pt5 = Pt42 * tau_lpt ** (gamma / (gamma - 1))

        # ----------------------
        # Nozzle primary flow
        # ----------------------
        Tt8 = Tt5
        Pt8 = Pt5
        M8 = 1
        M8_ad = np.sqrt((2 / (gamma - 1)) * ((Pt8 / Ps0) ** ((gamma - 1) / gamma)) - 1)
        if M8_ad < 1:
            M8 = M8_ad
        V8 = M8 * np.sqrt(gamma * R * Tt8)  # simplified

        # ----------------------
        # Nozzle secondary flow
        # ----------------------
        Tt18 = Tt13
        Pt18 = Pt13
        M18 = 1
        M18_ad = np.sqrt((2 / (gamma - 1)) * ((Pt18 / Ps0) ** ((gamma - 1) / gamma)) - 1)
        if M18_ad < 1:
            M18 = M18_ad
        V18 = M18 * np.sqrt(gamma * R * Tt18)

        # ----------------------
        # Thrust calculations
        # ----------------------
        V0 = M0 * np.sqrt(gamma * R * Ts0)
        T_fp = (1 + f) * dm_fp * V8
        T_fs = dm_fs * V18
        T = T_fp + T_fs - dm * V0
        ST = T / dm
        TSFC = dm_f / ST

        # ----------------------
        # Power and efficiency
        # ----------------------
        Wprop = T * V0
        Wcin = ((1 + f) * dm_fp * V8 ** 2) / 2 + (dm_fs * V18 ** 2) / 2 - (dm * V0 ** 2) / 2
        Wth = dm_f * PCI

        nprop = Wprop / Wcin
        nth = Wcin / Wth
        nglob = Wprop / Wth

        # ----------------------
        # Section areas (simplified)
        # ----------------------
        Dr1 = np.sqrt(gamma / R) * M1 * (1 + ((gamma - 1) / 2) * M1 ** 2) ** (-(gamma + 1) / (2 * (gamma - 1)))
        Dr2 = np.sqrt(gamma / R) * M2 * (1 + ((gamma - 1) / 2) * M2 ** 2) ** (-(gamma + 1) / (2 * (gamma - 1)))
        A1 = (dm * np.sqrt(Tt2)) / (Pt2 * Dr1)
        A2 = (dm * np.sqrt(Tt2)) / (Pt2 * Dr2)
        A8 = ((1 + f) * dm_fp * np.sqrt(Tt8)) / Pt8
        A18 = (dm_fs * np.sqrt(Tt18)) / Pt18

        # ----------------------
        # Save results
        # ----------------------
        results["T_list"].append(T)
        results["T_fp_list"].append(T_fp)
        results["T_fs_list"].append(T_fs)
        results["ST_list"].append(ST)
        results["TSFC_list"].append(TSFC)
        results["Wprop_list"].append(Wprop)
        results["Wcin_list"].append(Wcin)
        results["Wth_list"].append(Wth)
        results["nprop_list"].append(nprop)
        results["nth_list"].append(nth)
        results["nglob_list"].append(nglob)
        results["A1_list"].append(A1)
        results["A2_list"].append(A2)
        results["A8_list"].append(A8)
        results["A18_list"].append(A18)

    return results
