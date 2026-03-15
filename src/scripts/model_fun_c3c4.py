from types import SimpleNamespace

import numpy as np


def _field_any(v, *names):
    if isinstance(v, dict):
        for name in names:
            if name in v:
                return v[name]
    else:
        for name in names:
            if hasattr(v, name):
                return getattr(v, name)

    raise AttributeError(f"Missing required field. Tried: {names!r}")


def _normalize_pathway_option(pathway_option):
    if isinstance(pathway_option, (list, tuple)) and len(pathway_option) == 1:
        return pathway_option[0]
    return pathway_option


def _all_true(value):
    value = np.asarray(value)
    return value.size != 0 and bool(np.all(value))


_MODEL_OUTPUT_EXCLUDE = {
    "v",
    "Absorptance",
    "Cytbf_density",
    "Rubisco_density",
    "Resp_scalar",
    "k_F",
    "k_D",
    "k_P1",
    "k_N1",
    "k_P2",
    "k_U2",
    "k_q",
    "ATP_e_ratio_linear",
    "ATP_e_ratio_cyclic",
    "k_cat_CO2",
    "k_cat_O2",
    "K_m_CO2",
    "K_m_O2",
    "solve_C3C4_cc",
    "solve_C3C4_cj",
    "solve_C3C4_jc",
    "solve_C3C4_jj",
    "solve_C4_cc",
    "solve_C4_cj",
    "solve_C4_jc",
    "solve_C4_jj",
    "epsilon_PSI",
    "epsilon_PSII",
}


def build_model_output(values, exclude=None):
    exclude_names = set(_MODEL_OUTPUT_EXCLUDE)
    if exclude is not None:
        exclude_names.update(exclude)

    workspace = {
        name: value
        for name, value in values.items()
        if name not in exclude_names
    }
    return SimpleNamespace(**workspace)


def compute_fluorescence_outputs(state):
    state = dict(state)

    Cytbf_density_m = state["Cytbf_density"] * (1 - state["Cytbf_fraction_s"])
    Rubisco_density_m = state["Rubisco_density"] * (1 - state["Rubisco_fraction_s"])
    Cytbf_active_m_actual = state["J_PSI_m_j"] / state["k_q"]
    phi_P1_m_actual = state["J_PSI_m_actual"] / (state["PPFD"] * state["a_PSI_m"])
    q_P1_m_actual = phi_P1_m_actual * ((state["k_P1"] + state["k_D"] + state["k_F"]) / state["k_P1"])
    phi_P2_m_actual = state["J_PSII_m_actual"] / (state["PPFD"] * state["a_PSII_m"])
    q_P2_m_actual = 1 - Cytbf_active_m_actual / Cytbf_density_m

    k_N2_m_actual = (
        (
            state["k_P2"] ** 2 * phi_P2_m_actual**2
            - 2 * state["k_P2"] ** 2 * phi_P2_m_actual * q_P2_m_actual
            + state["k_P2"] ** 2 * q_P2_m_actual**2
            - 4 * state["k_P2"] * state["k_U2"] * phi_P2_m_actual**2 * q_P2_m_actual
            + 2 * state["k_P2"] * state["k_U2"] * phi_P2_m_actual**2
            + 2 * state["k_P2"] * state["k_U2"] * phi_P2_m_actual * q_P2_m_actual
            + state["k_U2"] ** 2 * phi_P2_m_actual**2
        ) ** (1 / 2)
        - state["k_P2"] * phi_P2_m_actual
        + state["k_U2"] * phi_P2_m_actual
        + state["k_P2"] * q_P2_m_actual
    ) / (2 * phi_P2_m_actual) - state["k_F"] - state["k_U2"] - state["k_D"]

    if _all_true(state["V_p_max_m"] == 0):
        Cytbf_density_s = 0
        Rubisco_density_s = 0
        Cytbf_active_s_actual = 0
        phi_P1_s_actual = 0
        q_P1_s_actual = 0
        phi_P2_s_actual = 0
        q_P2_s_actual = 0
        k_N2_s_actual = 0

    if _all_true(state["V_c_max_s"] > 0):
        Cytbf_density_s = state["Cytbf_density"] * state["Cytbf_fraction_s"]
        Rubisco_density_s = state["Rubisco_density"] * state["Rubisco_fraction_s"]
        Cytbf_active_s_actual = state["J_PSI_s_jj"] / state["k_q"]
        phi_P1_s_actual = state["J_PSI_s_actual"] / (state["PPFD"] * state["a_PSI_s"])
        q_P1_s_actual = phi_P1_s_actual * ((state["k_P1"] + state["k_D"] + state["k_F"]) / state["k_P1"])
        phi_P2_s_actual = state["J_PSII_s_actual"] / (state["PPFD"] * state["a_PSII_s"])
        q_P2_s_actual = 1 - Cytbf_active_s_actual / Cytbf_density_s

        k_N2_s_actual = (
            (
                state["k_P2"] ** 2 * phi_P2_s_actual**2
                - 2 * state["k_P2"] ** 2 * phi_P2_s_actual * q_P2_s_actual
                + state["k_P2"] ** 2 * q_P2_s_actual**2
                - 4 * state["k_P2"] * state["k_U2"] * phi_P2_s_actual**2 * q_P2_s_actual
                + 2 * state["k_P2"] * state["k_U2"] * phi_P2_s_actual**2
                + 2 * state["k_P2"] * state["k_U2"] * phi_P2_s_actual * q_P2_s_actual
                + state["k_U2"] ** 2 * phi_P2_s_actual**2
            ) ** (1 / 2)
            - state["k_P2"] * phi_P2_s_actual
            + state["k_U2"] * phi_P2_s_actual
            + state["k_P2"] * q_P2_s_actual
        ) / (2 * phi_P2_s_actual) - state["k_F"] - state["k_U2"] - state["k_D"]

    phi_p2_m_actual = q_P2_m_actual * state["k_P2"] / (state["k_P2"] + k_N2_m_actual + state["k_D"] + state["k_F"] + state["k_U2"])
    phi_n2_m_actual = q_P2_m_actual * k_N2_m_actual / (state["k_P2"] + k_N2_m_actual + state["k_D"] + state["k_F"] + state["k_U2"]) + (1 - q_P2_m_actual) * k_N2_m_actual / (k_N2_m_actual + state["k_D"] + state["k_F"] + state["k_U2"])
    phi_d2_m_actual = q_P2_m_actual * state["k_D"] / (state["k_P2"] + k_N2_m_actual + state["k_D"] + state["k_F"] + state["k_U2"]) + (1 - q_P2_m_actual) * state["k_D"] / (k_N2_m_actual + state["k_D"] + state["k_F"] + state["k_U2"])
    phi_f2_m_actual = q_P2_m_actual * state["k_F"] / (state["k_P2"] + k_N2_m_actual + state["k_D"] + state["k_F"] + state["k_U2"]) + (1 - q_P2_m_actual) * state["k_F"] / (k_N2_m_actual + state["k_D"] + state["k_F"] + state["k_U2"])
    phi_u2_m_actual = q_P2_m_actual * state["k_U2"] / (state["k_P2"] + k_N2_m_actual + state["k_D"] + state["k_F"] + state["k_U2"]) + (1 - q_P2_m_actual) * state["k_U2"] / (k_N2_m_actual + state["k_D"] + state["k_F"] + state["k_U2"])
    phi_P2_m_actual = phi_p2_m_actual / (1 - phi_u2_m_actual)
    phi_N2_m_actual = phi_n2_m_actual / (1 - phi_u2_m_actual)
    phi_D2_m_actual = phi_d2_m_actual / (1 - phi_u2_m_actual)
    phi_F2_m_actual = phi_f2_m_actual / (1 - phi_u2_m_actual)

    phi_P1_m_actual = q_P1_m_actual * state["k_P1"] / (state["k_P1"] + state["k_D"] + state["k_F"])
    phi_N1_m_actual = (1 - q_P1_m_actual) * state["k_N1"] / (state["k_N1"] + state["k_D"] + state["k_F"])
    phi_D1_m_actual = q_P1_m_actual * state["k_D"] / (state["k_P1"] + state["k_D"] + state["k_F"]) + (1 - q_P1_m_actual) * state["k_D"] / (state["k_N1"] + state["k_D"] + state["k_F"])
    phi_F1_m_actual = q_P1_m_actual * state["k_F"] / (state["k_P1"] + state["k_D"] + state["k_F"]) + (1 - q_P1_m_actual) * state["k_F"] / (state["k_N1"] + state["k_D"] + state["k_F"])

    phi_p2_s_actual = q_P2_s_actual * state["k_P2"] / (state["k_P2"] + k_N2_s_actual + state["k_D"] + state["k_F"] + state["k_U2"])
    phi_n2_s_actual = q_P2_s_actual * k_N2_s_actual / (state["k_P2"] + k_N2_s_actual + state["k_D"] + state["k_F"] + state["k_U2"]) + (1 - q_P2_s_actual) * k_N2_s_actual / (k_N2_s_actual + state["k_D"] + state["k_F"] + state["k_U2"])
    phi_d2_s_actual = q_P2_s_actual * state["k_D"] / (state["k_P2"] + k_N2_s_actual + state["k_D"] + state["k_F"] + state["k_U2"]) + (1 - q_P2_s_actual) * state["k_D"] / (k_N2_s_actual + state["k_D"] + state["k_F"] + state["k_U2"])
    phi_f2_s_actual = q_P2_s_actual * state["k_F"] / (state["k_P2"] + k_N2_s_actual + state["k_D"] + state["k_F"] + state["k_U2"]) + (1 - q_P2_s_actual) * state["k_F"] / (k_N2_s_actual + state["k_D"] + state["k_F"] + state["k_U2"])
    phi_u2_s_actual = q_P2_s_actual * state["k_U2"] / (state["k_P2"] + k_N2_s_actual + state["k_D"] + state["k_F"] + state["k_U2"]) + (1 - q_P2_s_actual) * state["k_U2"] / (k_N2_s_actual + state["k_D"] + state["k_F"] + state["k_U2"])
    phi_P2_s_actual = phi_p2_s_actual / (1 - phi_u2_s_actual)
    phi_N2_s_actual = phi_n2_s_actual / (1 - phi_u2_s_actual)
    phi_D2_s_actual = phi_d2_s_actual / (1 - phi_u2_s_actual)
    phi_F2_s_actual = phi_f2_s_actual / (1 - phi_u2_s_actual)

    phi_P1_s_actual = q_P1_s_actual * state["k_P1"] / (state["k_P1"] + state["k_D"] + state["k_F"])
    phi_N1_s_actual = (1 - q_P1_s_actual) * state["k_N1"] / (state["k_N1"] + state["k_D"] + state["k_F"])
    phi_D1_s_actual = q_P1_s_actual * state["k_D"] / (state["k_P1"] + state["k_D"] + state["k_F"]) + (1 - q_P1_s_actual) * state["k_D"] / (state["k_N1"] + state["k_D"] + state["k_F"])
    phi_F1_s_actual = q_P1_s_actual * state["k_F"] / (state["k_P1"] + state["k_D"] + state["k_F"]) + (1 - q_P1_s_actual) * state["k_F"] / (state["k_N1"] + state["k_D"] + state["k_F"])

    F_s_m_actual = state["a_PSII_m"] * phi_F2_m_actual * state["epsilon_PSII"] + state["a_PSI_m"] * phi_F1_m_actual * state["epsilon_PSI"]
    F_m_m_actual = state["a_PSII_m"] * state["k_F"] / (state["k_D"] + state["k_F"]) * state["epsilon_PSII"] + state["a_PSI_m"] * state["k_F"] / (state["k_N1"] + state["k_D"] + state["k_F"]) * state["epsilon_PSI"]
    F_o_m_actual = state["a_PSII_m"] * state["k_F"] / (state["k_P2"] + state["k_D"] + state["k_F"]) * state["epsilon_PSII"] + state["a_PSI_m"] * state["k_F"] / (state["k_P1"] + state["k_D"] + state["k_F"]) * state["epsilon_PSI"]
    F_m_prime_m_actual = state["a_PSII_m"] * state["k_F"] / (k_N2_m_actual + state["k_D"] + state["k_F"]) * state["epsilon_PSII"] + state["a_PSI_m"] * state["k_F"] / (state["k_N1"] + state["k_D"] + state["k_F"]) * state["epsilon_PSI"]
    F_o_prime_m_actual = state["a_PSII_m"] * state["k_F"] / (state["k_P2"] + k_N2_m_actual + state["k_D"] + state["k_F"]) * state["epsilon_PSII"] + state["a_PSI_m"] * state["k_F"] / (state["k_P1"] + state["k_D"] + state["k_F"]) * state["epsilon_PSI"]

    F_s_s_actual = state["a_PSII_s"] * phi_F2_s_actual * state["epsilon_PSII"] + state["a_PSI_s"] * phi_F1_s_actual * state["epsilon_PSI"]
    F_m_s_actual = state["a_PSII_s"] * state["k_F"] / (state["k_D"] + state["k_F"]) * state["epsilon_PSII"] + state["a_PSI_s"] * state["k_F"] / (state["k_N1"] + state["k_D"] + state["k_F"]) * state["epsilon_PSI"]
    F_o_s_actual = state["a_PSII_s"] * state["k_F"] / (state["k_P2"] + state["k_D"] + state["k_F"]) * state["epsilon_PSII"] + state["a_PSI_s"] * state["k_F"] / (state["k_P1"] + state["k_D"] + state["k_F"]) * state["epsilon_PSI"]
    F_m_prime_s_actual = state["a_PSII_s"] * state["k_F"] / (k_N2_s_actual + state["k_D"] + state["k_F"]) * state["epsilon_PSII"] + state["a_PSI_s"] * state["k_F"] / (state["k_N1"] + state["k_D"] + state["k_F"]) * state["epsilon_PSI"]
    F_o_prime_s_actual = state["a_PSII_s"] * state["k_F"] / (state["k_P2"] + k_N2_s_actual + state["k_D"] + state["k_F"]) * state["epsilon_PSII"] + state["a_PSI_s"] * state["k_F"] / (state["k_P1"] + state["k_D"] + state["k_F"]) * state["epsilon_PSI"]

    F_s = F_s_m_actual + F_s_s_actual
    F_m = F_m_m_actual + F_m_s_actual
    F_o = F_o_m_actual + F_o_s_actual
    F_m_prime = F_m_prime_m_actual + F_m_prime_s_actual
    F_o_prime = F_o_prime_m_actual + F_o_prime_s_actual

    qP = (F_m_prime - F_s) / (F_m_prime - F_o_prime)
    qL = (F_m_prime - F_s) * F_o_prime / ((F_m_prime - F_o_prime) * F_s)
    NPQ = F_m / F_m_prime - 1
    Phi_PSII = 1 - F_s / F_m_prime
    Phi_NPQ = F_s * (1 / F_m_prime - 1 / F_m)
    Phi_f_d = F_s / F_m

    return {
        "Cytbf_density_m": Cytbf_density_m,
        "Rubisco_density_m": Rubisco_density_m,
        "Cytbf_active_m_actual": Cytbf_active_m_actual,
        "phi_P1_m_actual": phi_P1_m_actual,
        "q_P1_m_actual": q_P1_m_actual,
        "q_P2_m_actual": q_P2_m_actual,
        "k_N2_m_actual": k_N2_m_actual,
        "Cytbf_density_s": Cytbf_density_s,
        "Rubisco_density_s": Rubisco_density_s,
        "Cytbf_active_s_actual": Cytbf_active_s_actual,
        "phi_P1_s_actual": phi_P1_s_actual,
        "q_P1_s_actual": q_P1_s_actual,
        "phi_P2_s_actual": phi_P2_s_actual,
        "q_P2_s_actual": q_P2_s_actual,
        "k_N2_s_actual": k_N2_s_actual,
        "phi_p2_m_actual": phi_p2_m_actual,
        "phi_n2_m_actual": phi_n2_m_actual,
        "phi_d2_m_actual": phi_d2_m_actual,
        "phi_f2_m_actual": phi_f2_m_actual,
        "phi_u2_m_actual": phi_u2_m_actual,
        "phi_P2_m_actual": phi_P2_m_actual,
        "phi_N2_m_actual": phi_N2_m_actual,
        "phi_D2_m_actual": phi_D2_m_actual,
        "phi_F2_m_actual": phi_F2_m_actual,
        "phi_N1_m_actual": phi_N1_m_actual,
        "phi_D1_m_actual": phi_D1_m_actual,
        "phi_F1_m_actual": phi_F1_m_actual,
        "phi_p2_s_actual": phi_p2_s_actual,
        "phi_n2_s_actual": phi_n2_s_actual,
        "phi_d2_s_actual": phi_d2_s_actual,
        "phi_f2_s_actual": phi_f2_s_actual,
        "phi_u2_s_actual": phi_u2_s_actual,
        "phi_P2_s_actual": phi_P2_s_actual,
        "phi_N2_s_actual": phi_N2_s_actual,
        "phi_D2_s_actual": phi_D2_s_actual,
        "phi_F2_s_actual": phi_F2_s_actual,
        "phi_N1_s_actual": phi_N1_s_actual,
        "phi_D1_s_actual": phi_D1_s_actual,
        "phi_F1_s_actual": phi_F1_s_actual,
        "F_s_m_actual": F_s_m_actual,
        "F_m_m_actual": F_m_m_actual,
        "F_o_m_actual": F_o_m_actual,
        "F_m_prime_m_actual": F_m_prime_m_actual,
        "F_o_prime_m_actual": F_o_prime_m_actual,
        "F_s_s_actual": F_s_s_actual,
        "F_m_s_actual": F_m_s_actual,
        "F_o_s_actual": F_o_s_actual,
        "F_m_prime_s_actual": F_m_prime_s_actual,
        "F_o_prime_s_actual": F_o_prime_s_actual,
        "F_s": F_s,
        "F_m": F_m,
        "F_o": F_o,
        "F_m_prime": F_m_prime,
        "F_o_prime": F_o_prime,
        "qP": qP,
        "qL": qL,
        "NPQ": NPQ,
        "Phi_PSII": Phi_PSII,
        "Phi_NPQ": Phi_NPQ,
        "Phi_f_d": Phi_f_d,
    }


def model_fun_c3c4(v):
    from .model_fun_c3 import model_fun_c3
    from .model_fun_nadp_me_c4 import model_fun_nadp_me_c4
    from .model_fun_type_i_c3c4 import model_fun_type_i_c3c4

    pathway_option = _normalize_pathway_option(_field_any(v, "pathway_option", "pathway_opt"))
    pathway_dispatch = {
        "C3": model_fun_c3,
        "Type-I-C3-C4": model_fun_type_i_c3c4,
        "NADP-ME-C4": model_fun_nadp_me_c4,
    }

    try:
        pathway_model_fun = pathway_dispatch[pathway_option]
    except KeyError as exc:
        supported = ", ".join(sorted(pathway_dispatch))
        raise ValueError(f"Unsupported pathway_option {pathway_option!r}. Supported pathways: {supported}") from exc

    return pathway_model_fun(v)
