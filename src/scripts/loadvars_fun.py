def _field(v, name):
    if isinstance(v, dict):
        return v[name]
    return getattr(v, name)


def loadvars_fun(v):
    pathway_opt = _field(v, "pathway_opt")
    Q = _field(v, "Q")
    T = _field(v, "T")
    P = _field(v, "P")
    O_m = _field(v, "O_m")
    C_m = _field(v, "C_m")
    Abs = _field(v, "Abs")
    abs_frac = _field(v, "abs_frac")
    a2_m_frac = _field(v, "a2_m_frac")
    a2_s_frac = _field(v, "a2_s_frac")
    vq_frac = _field(v, "vq_frac")
    vc_frac = _field(v, "vc_frac")
    CB6F = _field(v, "CB6F")
    RUB = _field(v, "RUB")
    Rdsc = _field(v, "Rdsc")
    Vpmax = _field(v, "Vpmax")
    gbs = _field(v, "gbs")
    gbso = _field(v, "gbso")
    Kf = _field(v, "Kf")
    Kd = _field(v, "Kd")
    Kp1 = _field(v, "Kp1")
    Kp2 = _field(v, "Kp2")
    Kn1 = _field(v, "Kn1")
    Ku2 = _field(v, "Ku2")
    kq = _field(v, "kq")
    nl = _field(v, "nl")
    nc = _field(v, "nc")
    kc = _field(v, "kc")
    ko = _field(v, "ko")
    Kc = _field(v, "Kc")
    Ko = _field(v, "Ko")
    Kp = _field(v, "Kp")
    eps1 = _field(v, "eps1")
    eps2 = _field(v, "eps2")
    c3c4_solve_cc = _field(v, "c3c4_solve_cc")
    c3c4_solve_cj = _field(v, "c3c4_solve_cj")
    c3c4_solve_jc = _field(v, "c3c4_solve_jc")
    c3c4_solve_jj = _field(v, "c3c4_solve_jj")
    c4_solve_cc = _field(v, "c4_solve_cc")
    c4_solve_cj = _field(v, "c4_solve_cj")
    c4_solve_jc = _field(v, "c4_solve_jc")
    c4_solve_jj = _field(v, "c4_solve_jj")

    return (
        pathway_opt,
        Q,
        T,
        P,
        O_m,
        C_m,
        Abs,
        abs_frac,
        a2_m_frac,
        a2_s_frac,
        eps1,
        eps2,
        CB6F,
        vq_frac,
        RUB,
        Rdsc,
        vc_frac,
        Vpmax,
        gbs,
        gbso,
        Kf,
        Kd,
        Kp1,
        Kn1,
        Kp2,
        Ku2,
        kq,
        nl,
        nc,
        kc,
        ko,
        Kc,
        Ko,
        Kp,
        c3c4_solve_cc,
        c3c4_solve_cj,
        c3c4_solve_jc,
        c3c4_solve_jj,
        c4_solve_cc,
        c4_solve_cj,
        c4_solve_jc,
        c4_solve_jj,
    )
