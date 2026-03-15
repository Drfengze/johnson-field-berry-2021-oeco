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


def loadvars_fun(v):
    # Legacy MATLAB-style aliases remain readable here so older inputs still load
    # while the rest of the Python code uses the standardized names.
    pathway_option = _field_any(v, "pathway_option", "pathway_opt")
    PPFD = _field_any(v, "PPFD", "Q")
    Temp = _field_any(v, "Temp", "T")
    Pressure = _field_any(v, "Pressure", "P")
    O2_m = _field_any(v, "O2_m", "O_m")
    CO2_m = _field_any(v, "CO2_m", "C_m")
    Absorptance = _field_any(v, "Absorptance", "Abs")
    Abs_fraction_s = _field_any(v, "Abs_fraction_s", "abs_frac")
    PSII_fraction = _field_any(v, "PSII_fraction", "a2_m_frac")
    PSII_fraction_s = _field_any(v, "PSII_fraction_s", "a2_s_frac")
    Cytbf_fraction_s = _field_any(v, "Cytbf_fraction_s", "vq_frac")
    Rubisco_fraction_s = _field_any(v, "Rubisco_fraction_s", "vc_frac")
    Cytbf_density = _field_any(v, "Cytbf_density", "CB6F")
    Rubisco_density = _field_any(v, "Rubisco_density", "RUB")
    Resp_scalar = _field_any(v, "Resp_scalar", "Rdsc")
    V_p_max = _field_any(v, "V_p_max", "Vpmax")
    g_bs_CO2 = _field_any(v, "g_bs_CO2", "gbs")
    g_bs_O2 = _field_any(v, "g_bs_O2", "gbso")
    k_F = _field_any(v, "k_F", "Kf")
    k_D = _field_any(v, "k_D", "Kd")
    k_P1 = _field_any(v, "k_P1", "Kp1")
    k_P2 = _field_any(v, "k_P2", "Kp2")
    k_N1 = _field_any(v, "k_N1", "Kn1")
    k_U2 = _field_any(v, "k_U2", "Ku2")
    k_q = _field_any(v, "k_q", "kq")
    ATP_e_ratio_linear = _field_any(v, "ATP_e_ratio_linear", "nl")
    ATP_e_ratio_cyclic = _field_any(v, "ATP_e_ratio_cyclic", "nc")
    k_cat_CO2 = _field_any(v, "k_cat_CO2", "kc")
    k_cat_O2 = _field_any(v, "k_cat_O2", "ko")
    K_m_CO2 = _field_any(v, "K_m_CO2", "Kc")
    K_m_O2 = _field_any(v, "K_m_O2", "Ko")
    K_m_PEPC_CO2 = _field_any(v, "K_m_PEPC_CO2", "Kp")
    epsilon_PSI = _field_any(v, "epsilon_PSI", "eps1")
    epsilon_PSII = _field_any(v, "epsilon_PSII", "eps2")
    solve_C3C4_cc = _field_any(v, "solve_C3C4_cc", "c3c4_solve_cc")
    solve_C3C4_cj = _field_any(v, "solve_C3C4_cj", "c3c4_solve_cj")
    solve_C3C4_jc = _field_any(v, "solve_C3C4_jc", "c3c4_solve_jc")
    solve_C3C4_jj = _field_any(v, "solve_C3C4_jj", "c3c4_solve_jj")
    solve_C4_cc = _field_any(v, "solve_C4_cc", "c4_solve_cc")
    solve_C4_cj = _field_any(v, "solve_C4_cj", "c4_solve_cj")
    solve_C4_jc = _field_any(v, "solve_C4_jc", "c4_solve_jc")
    solve_C4_jj = _field_any(v, "solve_C4_jj", "c4_solve_jj")

    return (
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
    )
