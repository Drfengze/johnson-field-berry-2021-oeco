import numpy as np

from .workspace2struct_fun import workspace2struct_fun


def _field(data, name):
    if isinstance(data, dict):
        return data[name]
    return getattr(data, name)


def configure_fun(data):
    pathway_option = "C3"

    PPFD = np.asarray(_field(data, "Qin"), dtype=float) / 1e6
    PPFD_u = "mol PAR m-2 s-1"
    PPFD_d = "PPFD incident on the leaf"

    Temp = np.asarray(_field(data, "Tin"), dtype=float)
    Temp_u = "degrees C"
    Temp_d = "Temp (leaf temperature)"

    Pressure = np.asarray(_field(data, "Pin"), dtype=float)
    Pressure_u = "bar"
    Pressure_d = "Pressure (total pressure)"

    CO2_m = np.asarray(_field(data, "Cin"), dtype=float) / 1e6
    CO2_m_u = "bar"
    CO2_m_d = "CO2_m (mesophyll CO2 partial pressure)"

    O2_m = np.asarray(_field(data, "Oin"), dtype=float)
    O2_m_u = "bar"
    O2_m_d = "O2_m (mesophyll O2 partial pressure)"

    Absorptance = 0.85
    Absorptance_u = "mol absorbed mol-1 incident"
    Absorptance_d = "Absorptance to PAR"

    Abs_fraction_s = 0.0
    Abs_fraction_s_u = "dimensionless"
    Abs_fraction_s_d = "Abs_fraction_s of total leaf absorptance to PAR"

    PSII_fraction = 0.52
    PSII_fraction_u = "dimensionless"
    PSII_fraction_d = "PSII_fraction of mesophyll absorptance to PAR"

    PSII_fraction_s = 0.0
    PSII_fraction_s_u = "dimensionless"
    PSII_fraction_s_d = "PSII_fraction_s of bundle sheath absorptance"

    Cytbf_fraction_s = 0.0
    Cytbf_fraction_s_u = "dimensionless"
    Cytbf_fraction_s_d = "Cytbf_fraction_s of total Cyt b6f density"

    Rubisco_fraction_s = 0.0
    Rubisco_fraction_s_u = "dimensionless"
    Rubisco_fraction_s_d = "Rubisco_fraction_s of total Rubisco density"

    Cytbf_density = (350 / 300) / 1e6
    Cytbf_density_u = "mol sites m-2"
    Cytbf_density_d = "Cytbf_density"

    Rubisco_density = (100 / 3.6) / 1e6
    Rubisco_density_u = "mol sites m-2"
    Rubisco_density_d = "Rubisco_density"

    Resp_scalar = 0.01
    Resp_scalar_u = "dimensionless"
    Resp_scalar_d = "Resp_scalar for mitochondrial dark respiration"

    V_p_max = 0.0
    V_p_max_u = "mol CO2 m-2 s-1"
    V_p_max_d = "V_p_max of PEP carboxylase"

    g_bs_CO2 = 0.003
    g_bs_CO2_u = "mol CO2 m-2 s-1 bar-1"
    g_bs_CO2_d = "g_bs_CO2 (bundle sheath conductance to CO2)"

    g_bs_O2 = g_bs_CO2 * 0.047
    g_bs_O2_u = "mol O2 m-2 s-1 bar-1"
    g_bs_O2_d = "g_bs_O2 (bundle sheath conductance to O2)"

    k_F = 0.05e09
    k_F_u = "s-1"
    k_F_d = "Rate constant for fluorescence at PSII and PSI"

    k_D = 0.55e09
    k_D_u = "s-1"
    k_D_d = "Rate constant for constitutive heat loss at PSII and PSI"

    k_P1 = 14.5e09
    k_P1_u = "s-1"
    k_P1_d = "Rate constant for photochemistry at PSI"

    k_N1 = 14.5e09
    k_N1_u = "s-1"
    k_N1_d = "Rate constant for regulated heat loss at PSI"

    k_P2 = 4.5e09
    k_P2_u = "s-1"
    k_P2_d = "Rate constant for photochemistry at PSII"

    k_U2 = 0e09
    k_U2_u = "s-1"
    k_U2_d = "Rate constant for exciton sharing at PSII"

    k_q = 300
    k_q_u = "mol e-1 mol sites-1 s-1"
    k_q_d = "k_q for Cyt b6f acting on PQH2"

    ATP_e_ratio_linear = 0.75
    ATP_e_ratio_linear_u = "ATP/e-"
    ATP_e_ratio_linear_d = "ATP_e_ratio_linear"

    ATP_e_ratio_cyclic = 1.00
    ATP_e_ratio_cyclic_u = "ATP/e-"
    ATP_e_ratio_cyclic_d = "ATP_e_ratio_cyclic"

    k_cat_CO2 = 3.6
    k_cat_CO2_u = "mol CO2 mol sites-1 s-1"
    k_cat_CO2_d = "k_cat_CO2 for Rubisco"

    k_cat_O2 = 3.6 * 0.27
    k_cat_O2_u = "mol O2 mol sites-1 s-1"
    k_cat_O2_d = "k_cat_O2 for Rubisco"

    K_m_CO2 = 260 / 1e6
    K_m_CO2_u = "bar"
    K_m_CO2_d = "K_m_CO2 for Rubisco"

    K_m_O2 = 179000 / 1e6
    K_m_O2_u = "bar"
    K_m_O2_d = "K_m_O2 for Rubisco"

    K_m_PEPC_CO2 = 80 / 1e6
    K_m_PEPC_CO2_u = "bar"
    K_m_PEPC_CO2_d = "K_m_PEPC_CO2 for PEPC"

    epsilon_PSI = 0
    epsilon_PSI_u = "mol PSI F to detector mol-1 PSI F emitted"
    epsilon_PSI_d = "epsilon_PSI transfer function"

    epsilon_PSII = 1
    epsilon_PSII_u = "mol PSII F to detector mol-1 PSII F emitted"
    epsilon_PSII_d = "epsilon_PSII transfer function"

    solve_C3C4_cc = "NAN"
    solve_C3C4_cj = "NAN"
    solve_C3C4_jc = "NAN"
    solve_C3C4_jj = "NAN"
    solve_C4_cc = "NAN"
    solve_C4_cj = "NAN"
    solve_C4_jc = "NAN"
    solve_C4_jj = "NAN"

    return workspace2struct_fun(exclude={"data"})
