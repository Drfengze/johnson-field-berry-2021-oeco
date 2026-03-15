import numpy as np

from .loadvars_fun import loadvars_fun
from .model_fun_c3c4 import _all_true, build_model_output, compute_fluorescence_outputs


def _strcmp(pathway_option, target):
    if isinstance(pathway_option, (list, tuple)) and len(pathway_option) == 1:
        pathway_option = pathway_option[0]
    return pathway_option == target


def _min_with_index(values):
    arrays = np.broadcast_arrays(*[np.asarray(value, dtype=float) for value in values])
    stacked = np.stack(arrays, axis=-1)
    return np.min(stacked, axis=-1), np.argmin(stacked, axis=-1) + 1


def _get_optional(v, name, default):
    if isinstance(v, dict):
        return v.get(name, default)
    return getattr(v, name, default)


def _smooth_transition(limit1, limit2, theta):
    limit1, limit2 = np.broadcast_arrays(np.asarray(limit1, dtype=float), np.asarray(limit2, dtype=float))
    discriminant = np.maximum((limit1 + limit2) ** 2 - 4 * theta * limit1 * limit2, 0)
    sqrt_term = np.sqrt(discriminant)
    root1 = ((limit1 + limit2) + sqrt_term) / (2 * theta)
    root2 = ((limit1 + limit2) - sqrt_term) / (2 * theta)
    return np.minimum(root1, root2)


def _resolve_mesophyll_cross_sections(
    absorptance,
    psii_fraction,
    cytbf_density,
    ppfd,
    eta,
    k_q,
    k_D,
    k_F,
    k_P1,
    k_P2,
    k_U2,
    alpha_option,
    solve_cross_sections,
):
    a_PSII = absorptance * psii_fraction
    a_PSI = absorptance - a_PSII

    if alpha_option != "dynamic" or solve_cross_sections is None:
        return a_PSII, a_PSI

    if np.all(np.asarray(absorptance, dtype=float) == 0) or np.all(np.asarray(cytbf_density, dtype=float) == 0):
        return a_PSII, a_PSI

    phi_P1_max = k_P1 / (k_P1 + k_D + k_F)
    a_PSII = solve_cross_sections(
        absorptance,
        cytbf_density,
        k_D,
        k_F,
        k_P2,
        k_U2,
        ppfd,
        eta,
        k_q,
        phi_P1_max,
    )
    a_PSI = absorptance - a_PSII
    return a_PSII, a_PSI


def model_fun_c3(v):
    (
        pathway_option,
        PPFD,
        Temp,
        Pressure,
        O2_m,
        CO2_m,
        Absorptance,
        Abs_fraction_s,
        PSII_fraction,
        PSII_fraction_s,
        epsilon_PSI,
        epsilon_PSII,
        Cytbf_density,
        Cytbf_fraction_s,
        Rubisco_density,
        Resp_scalar,
        Rubisco_fraction_s,
        V_p_max,
        g_bs_CO2,
        g_bs_O2,
        k_F,
        k_D,
        k_P1,
        k_N1,
        k_P2,
        k_U2,
        k_q,
        ATP_e_ratio_linear,
        ATP_e_ratio_cyclic,
        k_cat_CO2,
        k_cat_O2,
        K_m_CO2,
        K_m_O2,
        K_m_PEPC_CO2,
        solve_C3C4_cc,
        solve_C3C4_cj,
        solve_C3C4_jc,
        solve_C3C4_jj,
        solve_C4_cc,
        solve_C4_cj,
        solve_C4_jc,
        solve_C4_jj,
    ) = loadvars_fun(v)

    PPFD = np.asarray(PPFD, dtype=float)
    Temp = np.asarray(Temp, dtype=float)
    Pressure = np.asarray(Pressure, dtype=float)
    O2_m = np.asarray(O2_m, dtype=float)
    CO2_m = np.asarray(CO2_m, dtype=float)

    alpha_option = _get_optional(v, "alpha_option", "static")
    solve_cross_sections = _get_optional(v, "solve_cross_sections", None)
    limitation_transition = _get_optional(v, "limitation_transition", "hard")
    theta_curvature = float(_get_optional(v, "theta_curvature", 0.95))
    if limitation_transition not in {"hard", "smooth"}:
        raise ValueError(
            "model_fun_c3 expects limitation_transition to be either 'hard' or 'smooth'."
        )

    Absorptance_m = Absorptance * (1 - Abs_fraction_s)
    Absorptance_s = Absorptance * Abs_fraction_s
    a_PSII_s = Absorptance_s * PSII_fraction_s
    a_PSI_s = Absorptance_s - a_PSII_s
    Cytbf_density_m = Cytbf_density * (1 - Cytbf_fraction_s)

    Temp_K = Temp + 273.15
    Temp_ref_K = 25 + 273.15
    R = 0.008314

    Ha = 37
    k_q = k_q * np.exp(Ha / R * (1 / Temp_ref_K - 1 / Temp_K))
    V_q_max = Cytbf_density * k_q
    V_q_max_m = V_q_max * (1 - Cytbf_fraction_s)
    V_q_max_s = V_q_max * Cytbf_fraction_s

    En = 0.710
    Hd = 220
    ATP_e_ratio_linear = ATP_e_ratio_linear * (1 + np.exp((Temp_ref_K * En - Hd) / (R * Temp_ref_K))) / (
        1 + np.exp((Temp_K * En - Hd) / (R * Temp_K))
    )
    ATP_e_ratio_cyclic = ATP_e_ratio_cyclic * (1 + np.exp((Temp_ref_K * En - Hd) / (R * Temp_ref_K))) / (
        1 + np.exp((Temp_K * En - Hd) / (R * Temp_K))
    )

    R_d = Rubisco_density * k_cat_CO2 * Resp_scalar
    Ha = 66
    R_d = R_d * np.exp(Ha / R * (1 / Temp_ref_K - 1 / Temp_K))
    R_d_m = R_d * (1 - Abs_fraction_s)
    R_d_s = R_d * Abs_fraction_s

    Specificity = (k_cat_CO2 / K_m_CO2) * (K_m_O2 / k_cat_O2)
    Ha = 23
    Specificity = 1 / (1 / Specificity * (np.exp(Ha / R * (1 / Temp_ref_K - 1 / Temp_K))))
    eta = (
        1
        - (ATP_e_ratio_linear / ATP_e_ratio_cyclic)
        + (3 + 7 * O2_m / (2 * Specificity * CO2_m))
        / ((4 + 4 * O2_m / (Specificity * CO2_m)) * ATP_e_ratio_cyclic)
    )
    a_PSII_m, a_PSI_m = _resolve_mesophyll_cross_sections(
        Absorptance_m,
        PSII_fraction,
        Cytbf_density_m,
        PPFD,
        eta,
        k_q,
        k_D,
        k_F,
        k_P1,
        k_P2,
        k_U2,
        alpha_option,
        solve_cross_sections,
    )
    Ha = 59
    K_m_CO2 = K_m_CO2 * np.exp(Ha / R * (1 / Temp_ref_K - 1 / Temp_K))
    Ha = 36
    K_m_O2 = K_m_O2 * np.exp(Ha / R * (1 / Temp_ref_K - 1 / Temp_K))
    Ha = 58
    k_cat_CO2 = k_cat_CO2 * np.exp(Ha / R * (1 / Temp_ref_K - 1 / Temp_K))
    V_c_max = Rubisco_density * k_cat_CO2
    V_c_max_m = V_c_max * (1 - Rubisco_fraction_s)
    V_c_max_s = V_c_max * Rubisco_fraction_s

    Ha = 58
    V_p_max = V_p_max * np.exp(Ha / R * (1 / Temp_ref_K - 1 / Temp_K))
    V_p_max_m = V_p_max

    if not _strcmp(pathway_option, "C3"):
        raise ValueError(f"model_fun_c3 only supports pathway_option='C3', got {pathway_option!r}")

    if not _all_true(V_p_max_m == 0):
        raise ValueError("model_fun_c3 expects V_p_max to be zero for the C3 pathway.")

    if np.mean(np.asarray(V_c_max_s, dtype=float)) > 0:
        raise ValueError("model_fun_c3 expects zero bundle-sheath Rubisco capacity for the C3 pathway.")

    J_PSI_m_j = PPFD * V_q_max_m / (PPFD + V_q_max_m / (a_PSI_m * (k_P1 / (k_P1 + k_D + k_F))))
    J_PSII_m_j = J_PSI_m_j / eta
    V_c_m_j = J_PSII_m_j / (4 * (1 + O2_m / (Specificity * CO2_m)))
    V_o_m_j = V_c_m_j * O2_m / (Specificity * CO2_m)
    A_gross_m_j = V_c_m_j - V_o_m_j / 2
    A_net_m_j = A_gross_m_j - R_d_m
    V_p_m_j = 0

    V_c_m_c = CO2_m * V_c_max_m / (CO2_m + K_m_CO2 * (1 + O2_m / K_m_O2))
    V_o_m_c = V_c_m_c * O2_m / (Specificity * CO2_m)
    A_gross_m_c = V_c_m_c - V_o_m_c / 2
    A_net_m_c = A_gross_m_c - R_d_m
    J_PSII_m_c = A_gross_m_c * 4 * (1 + O2_m / (Specificity * CO2_m)) / (1 - O2_m / (2 * Specificity * CO2_m))
    J_PSI_m_c = J_PSII_m_c * eta
    V_p_m_c = 0
    V_g_m_j = 0
    V_g_m_c = 0

    J_PSI_s_jj = 0
    J_PSII_s_jj = 0
    J_PSI_s_jc = 0
    J_PSII_s_jc = 0
    J_PSI_s_cc = 0
    J_PSII_s_cc = 0
    J_PSI_s_cj = 0
    J_PSII_s_cj = 0
    A_net_s_jj = 0
    A_net_s_jc = 0
    A_net_s_cj = 0
    A_net_s_cc = 0
    A_gross_s_jj = 0
    A_gross_s_jc = 0
    A_gross_s_cj = 0
    A_gross_s_cc = 0
    CO2_s_jj = 0
    CO2_s_jc = 0
    CO2_s_cj = 0
    CO2_s_cc = 0
    O2_s_jj = 0
    O2_s_jc = 0
    O2_s_cj = 0
    O2_s_cc = 0
    V_c_s_jj = 0
    V_c_s_jc = 0
    V_c_s_cj = 0
    V_c_s_cc = 0
    V_o_s_jj = 0
    V_o_s_jc = 0
    V_o_s_cj = 0
    V_o_s_cc = 0
    J_PSI_s_j = 0
    J_PSII_s_j = 0
    J_PSI_s_c = 0
    J_PSII_s_c = 0
    which_J_PSI_s_j = 0
    which_J_PSI_s_c = 0
    CO2_s_j = 0
    CO2_s_c = 0
    O2_s_j = 0
    O2_s_c = 0
    A_net_s_j = 0
    A_net_s_c = 0
    J_PSI_s_actual = 0
    J_PSII_s_actual = 0
    A_gross_s_actual = 0
    A_net_s_actual = 0
    CO2_s_actual = 0
    O2_s_actual = 0
    Leak_CO2_s_actual = 0

    V_g_m_actual = 0
    V_p_m_actual = 0

    if limitation_transition == "smooth":
        J_PSI_m_actual = _smooth_transition(J_PSI_m_j, J_PSI_m_c, theta_curvature)
        J_PSII_m_actual = _smooth_transition(J_PSII_m_j, J_PSII_m_c, theta_curvature)
        A_net_m_actual = _smooth_transition(A_net_m_j, A_net_m_c, theta_curvature)
        which_J_PSI_m = np.zeros_like(J_PSI_m_actual, dtype=int)
    else:
        J_PSI_m_actual, which_J_PSI_m = _min_with_index([J_PSI_m_j, J_PSI_m_c])
        J_PSII_m_actual = J_PSII_m_j * (which_J_PSI_m == 1) + J_PSII_m_c * (which_J_PSI_m == 2)
        A_net_m_actual = A_net_m_j * (which_J_PSI_m == 1) + A_net_m_c * (which_J_PSI_m == 2)

    A_gross_m_actual = A_net_m_actual + R_d_m

    if _all_true(V_c_max_s > 0):
        J_PSI_s_j, which_J_PSI_s_j = _min_with_index([J_PSI_s_jj, J_PSI_s_jc])
        J_PSII_s_j = J_PSII_s_jj * (which_J_PSI_s_j == 1) + J_PSII_s_jc * (which_J_PSI_s_j == 2)
        CO2_s_j = CO2_s_jj * (which_J_PSI_s_j == 1) + CO2_s_jc * (which_J_PSI_s_j == 2)
        O2_s_j = O2_s_jj * (which_J_PSI_s_j == 1) + O2_s_jc * (which_J_PSI_s_j == 2)
        A_net_s_j = A_net_s_jj * (which_J_PSI_s_j == 1) + A_net_s_jc * (which_J_PSI_s_j == 2)

        J_PSI_s_c, which_J_PSI_s_c = _min_with_index([J_PSI_s_cj, J_PSI_s_cc])
        J_PSII_s_c = J_PSII_s_cj * (which_J_PSI_s_c == 1) + J_PSII_s_cc * (which_J_PSI_s_c == 2)
        CO2_s_c = CO2_s_cj * (which_J_PSI_s_c == 1) + CO2_s_cc * (which_J_PSI_s_c == 2)
        O2_s_c = O2_s_cj * (which_J_PSI_s_c == 1) + O2_s_cc * (which_J_PSI_s_c == 2)
        A_net_s_c = A_net_s_cj * (which_J_PSI_s_c == 1) + A_net_s_cc * (which_J_PSI_s_c == 2)

        J_PSI_s_actual = J_PSI_s_j * (which_J_PSI_m == 1) + J_PSI_s_c * (which_J_PSI_m == 2)
        J_PSII_s_actual = J_PSII_s_j * (which_J_PSI_m == 1) + J_PSII_s_c * (which_J_PSI_m == 2)
        CO2_s_actual = CO2_s_j * (which_J_PSI_m == 1) + CO2_s_c * (which_J_PSI_m == 2)
        O2_s_actual = O2_s_j * (which_J_PSI_m == 1) + O2_s_c * (which_J_PSI_m == 2)
        A_net_s_actual = A_net_s_j * (which_J_PSI_m == 1) + A_net_s_c * (which_J_PSI_m == 2)
        Leak_CO2_s_actual = g_bs_CO2 / Pressure * (CO2_s_actual - CO2_m)
        A_gross_s_actual = A_net_s_actual + R_d_s

    J_PSI_actual = J_PSI_m_actual + J_PSI_s_actual
    J_PSII_actual = J_PSII_m_actual + J_PSII_s_actual
    A_net_actual = A_net_m_actual + A_net_s_actual
    A_gross_actual = A_net_actual + R_d

    fluorescence_inputs = {
        "PPFD": PPFD,
        "V_p_max_m": V_p_max_m,
        "V_c_max_s": V_c_max_s,
        "Cytbf_density": Cytbf_density,
        "Cytbf_fraction_s": Cytbf_fraction_s,
        "Rubisco_density": Rubisco_density,
        "Rubisco_fraction_s": Rubisco_fraction_s,
        "J_PSI_m_j": J_PSI_m_j,
        "J_PSI_m_actual": J_PSI_m_actual,
        "J_PSII_m_actual": J_PSII_m_actual,
        "J_PSI_s_jj": J_PSI_s_jj,
        "J_PSI_s_actual": J_PSI_s_actual,
        "J_PSII_s_actual": J_PSII_s_actual,
        "a_PSI_m": a_PSI_m,
        "a_PSII_m": a_PSII_m,
        "a_PSI_s": a_PSI_s,
        "a_PSII_s": a_PSII_s,
        "k_q": k_q,
        "k_P1": k_P1,
        "k_N1": k_N1,
        "k_P2": k_P2,
        "k_U2": k_U2,
        "k_D": k_D,
        "k_F": k_F,
        "epsilon_PSI": epsilon_PSI,
        "epsilon_PSII": epsilon_PSII,
    }
    fluorescence_outputs = compute_fluorescence_outputs(fluorescence_inputs)
    model_state = dict(locals())
    model_state.update(fluorescence_outputs)
    model_state.pop("fluorescence_inputs", None)
    model_state.pop("fluorescence_outputs", None)
    return build_model_output(
        model_state,
        exclude={
            "alpha_option",
            "solve_cross_sections",
            "limitation_transition",
            "theta_curvature",
            "eta",
        },
    )
