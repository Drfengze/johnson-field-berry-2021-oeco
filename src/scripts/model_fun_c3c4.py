import numpy as np

from .loadvars_fun import loadvars_fun
from .workspace2struct_fun import workspace2struct_fun


def _strcmp(pathway_opt, target):
    if isinstance(pathway_opt, (list, tuple)) and len(pathway_opt) == 1:
        pathway_opt = pathway_opt[0]
    return pathway_opt == target


def _matlab_truth(value):
    value = np.asarray(value)
    return value.size != 0 and bool(np.all(value))


def _min_with_index(values):
    arrays = np.broadcast_arrays(*[np.asarray(value, dtype=float) for value in values])
    stacked = np.stack(arrays, axis=-1)
    return np.min(stacked, axis=-1), np.argmin(stacked, axis=-1) + 1


def model_fun_c3c4(v):
    (
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
    ) = loadvars_fun(v)

    Q = np.asarray(Q, dtype=float)
    T = np.asarray(T, dtype=float)
    P = np.asarray(P, dtype=float)
    O_m = np.asarray(O_m, dtype=float)
    C_m = np.asarray(C_m, dtype=float)

    Abs_m = Abs * (1 - abs_frac)
    Abs_s = Abs * abs_frac
    a2_m = Abs_m * a2_m_frac
    a1_m = Abs_m - a2_m
    a2_s = Abs_s * a2_s_frac
    a1_s = Abs_s - a2_s

    T_K = T + 273.15
    Tref_K = 25 + 273.15
    R = 0.008314

    Ha = 37
    kq = kq * np.exp(Ha / R * (1 / Tref_K - 1 / T_K))
    Vqmax = CB6F * kq
    Vqmax_m = Vqmax * (1 - vq_frac)
    Vqmax_s = Vqmax * vq_frac

    En = 0.710
    Hd = 220
    nl = nl * (1 + np.exp((Tref_K * En - Hd) / (R * Tref_K))) / (
        1 + np.exp((T_K * En - Hd) / (R * T_K))
    )
    nc = nc * (1 + np.exp((Tref_K * En - Hd) / (R * Tref_K))) / (
        1 + np.exp((T_K * En - Hd) / (R * T_K))
    )

    Rd = RUB * kc * Rdsc
    Ha = 66
    Rd = Rd * np.exp(Ha / R * (1 / Tref_K - 1 / T_K))
    Rd_m = Rd * (1 - abs_frac)
    Rd_s = Rd * abs_frac

    S = (kc / Kc) * (Ko / ko)
    Ha = 23
    S = 1 / (1 / S * (np.exp(Ha / R * (1 / Tref_K - 1 / T_K))))
    Ha = 59
    Kc = Kc * np.exp(Ha / R * (1 / Tref_K - 1 / T_K))
    Ha = 36
    Ko = Ko * np.exp(Ha / R * (1 / Tref_K - 1 / T_K))
    Ha = 58
    kc = kc * np.exp(Ha / R * (1 / Tref_K - 1 / T_K))
    Vcmax = RUB * kc
    Vcmax_m = Vcmax * (1 - vc_frac)
    Vcmax_s = Vcmax * vc_frac

    Ha = 58
    Vpmax = Vpmax * np.exp(Ha / R * (1 / Tref_K - 1 / T_K))
    Vpmax_m = Vpmax

    if _matlab_truth(Vpmax_m == 0):
        if _strcmp(pathway_opt, "C3") or _strcmp(pathway_opt, "Type-I-C3-C4"):
            JP700_mj = Q * Vqmax_m / (Q + Vqmax_m / (a1_m * (Kp1 / (Kp1 + Kd + Kf))))
            JP680_mj = JP700_mj / (
                1
                - (nl / nc)
                + (3 + 7 * O_m / (2 * S * C_m))
                / ((4 + 4 * O_m / (S * C_m)) * nc)
            )
            Vc_mj = JP680_mj / (4 * (1 + O_m / (S * C_m)))
            Vo_mj = Vc_mj * O_m / (S * C_m)
            Ag_mj = Vc_mj - Vo_mj / 2
            An_mj = Ag_mj - Rd_m
            Vp_mj = 0

            Vc_mc = C_m * Vcmax_m / (C_m + Kc * (1 + O_m / Ko))
            Vo_mc = Vc_mc * O_m / (S * C_m)
            Ag_mc = Vc_mc - Vo_mc / 2
            An_mc = Ag_mc - Rd_m
            JP680_mc = Ag_mc * 4 * (1 + O_m / (S * C_m)) / (1 - O_m / (2 * S * C_m))
            JP700_mc = JP680_mc * (
                1
                - (nl / nc)
                + (3 + 7 * O_m / (2 * S * C_m))
                / ((4 + 4 * O_m / (S * C_m)) * nc)
            )
            Vp_mc = 0

            if _matlab_truth(Vcmax_s == 0):
                Vg_mj = 0
                Vg_mc = 0
            else:
                Vg_mj = Vo_mj / 2
                Vg_mc = Vo_mc / 2

    if _matlab_truth(Vpmax_m > 0):
        if _strcmp(pathway_opt, "NADP-ME-C4"):
            JP700_mj = Q * Vqmax_m / (Q + Vqmax_m / (a1_m * (Kp1 / (Kp1 + Kd + Kf))))
            JP680_mj = JP700_mj / (1 - (nl / nc) + 1 / nc)
            Vp_mj = JP680_mj / 2

            Vp_mc = Vpmax_m * C_m / (Kp + C_m)
            JP680_mc = Vp_mc * 2
            JP700_mc = JP680_mc * (1 - (nl / nc) + 1 / nc)

            Vc_mj = 0
            Vo_mj = 0
            Vg_mj = 0
            Ag_mj = 0
            An_mj = -Rd_m
            Vc_mc = 0
            Vo_mc = 0
            Vg_mc = 0
            Ag_mc = 0
            An_mc = -Rd_m

    if _matlab_truth(Vcmax_s == 0):
        JP700_sjj = 0
        JP680_sjj = 0
        JP700_sjc = 0
        JP680_sjc = 0
        JP700_scc = 0
        JP680_scc = 0
        JP700_scj = 0
        JP680_scj = 0
        JP700_sj = 0
        JP680_sj = 0
        JP700_sc = 0
        JP680_sc = 0
        JP700_sa = 0
        JP680_sa = 0
        which_JP700_sj = 0
        which_JP700_sc = 0
        Ag_sa = 0
        An_sa = 0
        C_sa = 0
        O_sa = 0
        L_C_sa = 0

    if np.mean(np.asarray(Vcmax_s, dtype=float)) > 0 and np.mean(np.asarray(Vpmax_m, dtype=float)) == 0:
        JP700_sjj = Q * Vqmax_s / (Q + Vqmax_s / (a1_s * (Kp1 / (Kp1 + Kd + Kf))))
        An_sjj = c3c4_solve_jj(C_m, JP680_mj, JP700_sjj, O_m, P, Rd_s, S, gbs, gbso, nc, nl)
        C_sjj = C_m + (
            (JP680_mj * O_m / (2 * S * C_m) / (4 * (1 + O_m / (S * C_m)))) - An_sjj
        ) * P / gbs
        O_sjj = O_m + An_sjj * P / gbso
        Ag_sjj = An_sjj + Rd_s
        JP680_sjj = Ag_sjj * 4 * (1 + O_sjj / (S * C_sjj)) / (1 - O_sjj / (2 * S * C_sjj))
        Vc_sjj = JP680_sjj / (4 * (1 + O_sjj / (S * C_sjj)))
        Vo_sjj = Vc_sjj * O_sjj / (S * C_sjj)

        An_sjc = c3c4_solve_jc(C_m, JP680_mj, Kc, Ko, O_m, P, Rd_s, S, Vcmax_s, gbs, gbso)
        C_sjc = C_m + (
            (JP680_mj * O_m / (2 * S * C_m) / (4 * (1 + O_m / (S * C_m)))) - An_sjc
        ) * P / gbs
        O_sjc = O_m + An_sjc * P / gbso
        Ag_sjc = An_sjc + Rd_s
        JP680_sjc = Ag_sjc * 4 * (1 + O_sjc / (S * C_sjc)) / (1 - O_sjc / (2 * S * C_sjc))
        JP700_sjc = JP680_sjc * (
            1
            - (nl / nc)
            + (3 + 7 * O_sjc / (2 * S * C_sjc)) / ((4 + 4 * O_sjc / (S * C_sjc)) * nc)
        )
        Vc_sjc = JP680_sjc / (4 * (1 + O_sjc / (S * C_sjc)))
        Vo_sjc = Vc_sjc * O_sjc / (S * C_sjc)

        JP700_scj = JP700_sjj
        An_scj = c3c4_solve_cj(C_m, JP700_scj, Kc, Ko, O_m, P, Rd_s, S, Vcmax_m, gbs, gbso, nc, nl)
        C_scj = C_m + (
            (Vcmax_m * O_m / (2 * S) / (C_m + Kc * (1 + O_m / Ko))) - An_scj
        ) * P / gbs
        O_scj = O_m + An_scj * P / gbso
        Ag_scj = An_scj + Rd_s
        JP680_scj = Ag_scj * 4 * (1 + O_scj / (S * C_scj)) / (1 - O_scj / (2 * S * C_scj))
        Vc_scj = JP680_scj / (4 * (1 + O_scj / (S * C_scj)))
        Vo_scj = Vc_scj * O_scj / (S * C_scj)

        An_scc = c3c4_solve_cc(C_m, Kc, Ko, O_m, P, Rd_s, S, Vcmax_m, Vcmax_s, gbs, gbso)
        C_scc = C_m + (
            (Vcmax_m * O_m / (2 * S) / (C_m + Kc * (1 + O_m / Ko))) - An_scc
        ) * P / gbs
        O_scc = O_m + An_scc * P / gbso
        Ag_scc = An_scc + Rd_s
        JP680_scc = Ag_scc * 4 * (1 + O_scc / (S * C_scc)) / (1 - O_scc / (2 * S * C_scc))
        JP700_scc = JP680_scc * (
            1
            - (nl / nc)
            + (3 + 7 * O_scc / (2 * S * C_scc)) / ((4 + 4 * O_scc / (S * C_scc)) * nc)
        )
        Vc_scc = JP680_scc / (4 * (1 + O_scc / (S * C_scc)))
        Vo_scc = Vc_scc * O_scc / (S * C_scc)

    if _matlab_truth(Vpmax_m > 0):
        if _strcmp(pathway_opt, "NADP-ME-C4"):
            JP700_sjj = Q * Vqmax_s / (Q + Vqmax_s / (a1_s * (Kp1 / (Kp1 + Kd + Kf))))
            An_sjj = c4_solve_jj(C_m, JP700_mj, JP700_sjj, O_m, P, Rd_s, S, gbs, gbso, nc, nl)
            C_sjj = C_m + (Vp_mj - An_sjj) * P / gbs
            O_sjj = O_m + (An_sjj - Vp_mj / 2) * P / gbso
            Ag_sjj = An_sjj + Rd_s
            JP680_sjj = Ag_sjj * 4 * (1 + O_sjj / (S * C_sjj)) / (1 - O_sjj / (2 * S * C_sjj)) - 2 * Vp_mj

            An_sjc = c4_solve_jc(C_m, JP700_mj, Kc, Ko, O_m, P, Rd_s, S, Vcmax_s, gbs, gbso, nc, nl)
            C_sjc = C_m + (Vp_mj - An_sjc) * P / gbs
            O_sjc = O_m + (An_sjc - Vp_mj / 2) * P / gbso
            Ag_sjc = An_sjc + Rd_s
            JP680_sjc = Ag_sjc * 4 * (1 + O_sjc / (S * C_sjc)) / (1 - O_sjc / (2 * S * C_sjc)) - 2 * Vp_mj
            JP700_sjc = JP680_sjc * (
                1
                - (nl / nc)
                + (3 + 7 * O_sjc / (2 * S * C_sjc)) / ((4 + 4 * O_sjc / (S * C_sjc)) * nc)
            ) + (3 + 7 * O_sjc / (2 * S * C_sjc)) / ((4 + 4 * O_sjc / (S * C_sjc)) * nc) * (Vp_mj / 2)

            JP700_scj = Q * Vqmax_s / (Q + Vqmax_s / (a1_s * (Kp1 / (Kp1 + Kd + Kf))))
            An_scj = c4_solve_cj(C_m, JP700_scj, Kp, O_m, P, Rd_s, S, Vpmax_m, gbs, gbso, nc, nl)
            C_scj = C_m + (Vp_mc - An_scj) * P / gbs
            O_scj = O_m + (An_scj - Vp_mc / 2) * P / gbso
            Ag_scj = An_scj + Rd_s
            JP680_scj = Ag_scj * 4 * (1 + O_scj / (S * C_scj)) / (1 - O_scj / (2 * S * C_scj)) - 2 * Vp_mc

            An_scc = c4_solve_cc(C_m, Kc, Ko, Kp, O_m, P, Rd_s, S, Vcmax_s, Vpmax_m, gbs, gbso)
            C_scc = C_m + (Vp_mc - An_scc) * P / gbs
            O_scc = O_m + (An_scc - Vp_mc / 2) * P / gbso
            Ag_scc = An_scc + Rd_s
            JP680_scc = Ag_scc * 4 * (1 + O_scc / (S * C_scc)) / (1 - O_scc / (2 * S * C_scc)) - 2 * Vp_mc
            JP700_scc = JP680_scc * (
                1
                - (nl / nc)
                + (3 + 7 * O_scc / (2 * S * C_scc)) / ((4 + 4 * O_scc / (S * C_scc)) * nc)
            ) + (3 + 7 * O_scc / (2 * S * C_scc)) / ((4 + 4 * O_scc / (S * C_scc)) * nc) * (Vp_mc / 2)

    JP700_ma, which_JP700_ma = _min_with_index([JP700_mj, JP700_mc])
    JP680_ma = JP680_mj * (which_JP700_ma == 1) + JP680_mc * (which_JP700_ma == 2)
    Vg_ma = Vg_mj * (which_JP700_ma == 1) + Vg_mc * (which_JP700_ma == 2)
    Vp_ma = Vp_mj * (which_JP700_ma == 1) + Vp_mc * (which_JP700_ma == 2)
    An_ma = An_mj * (which_JP700_ma == 1) + An_mc * (which_JP700_ma == 2)
    Ag_ma = An_ma + Rd_m

    if _matlab_truth(Vcmax_s > 0):
        JP700_sj, which_JP700_sj = _min_with_index([JP700_sjj, JP700_sjc])
        JP680_sj = JP680_sjj * (which_JP700_sj == 1) + JP680_sjc * (which_JP700_sj == 2)
        C_sj = C_sjj * (which_JP700_sj == 1) + C_sjc * (which_JP700_sj == 2)
        O_sj = O_sjj * (which_JP700_sj == 1) + O_sjc * (which_JP700_sj == 2)
        An_sj = An_sjj * (which_JP700_sj == 1) + An_sjc * (which_JP700_sj == 2)

        JP700_sc, which_JP700_sc = _min_with_index([JP700_scj, JP700_scc])
        JP680_sc = JP680_scj * (which_JP700_sc == 1) + JP680_scc * (which_JP700_sc == 2)
        C_sc = C_scj * (which_JP700_sc == 1) + C_scc * (which_JP700_sc == 2)
        O_sc = O_scj * (which_JP700_sc == 1) + O_scc * (which_JP700_sc == 2)
        An_sc = An_scj * (which_JP700_sc == 1) + An_scc * (which_JP700_sc == 2)

        JP700_sa = JP700_sj * (which_JP700_ma == 1) + JP700_sc * (which_JP700_ma == 2)
        JP680_sa = JP680_sj * (which_JP700_ma == 1) + JP680_sc * (which_JP700_ma == 2)
        C_sa = C_sj * (which_JP700_ma == 1) + C_sc * (which_JP700_ma == 2)
        O_sa = O_sj * (which_JP700_ma == 1) + O_sc * (which_JP700_ma == 2)
        An_sa = An_sj * (which_JP700_ma == 1) + An_sc * (which_JP700_ma == 2)
        L_C_sa = gbs / P * (C_sa - C_m)
        Ag_sa = An_sa + Rd_s

    JP700_a = JP700_ma + JP700_sa
    JP680_a = JP680_ma + JP680_sa
    An_a = An_ma + An_sa
    Ag_a = An_a + Rd

    CB6F_m = CB6F * (1 - vq_frac)
    RUB_m = RUB * (1 - vc_frac)
    CB6F_ma = JP700_mj / kq
    phi1P_ma = JP700_ma / (Q * a1_m)
    q1_ma = phi1P_ma * ((Kp1 + Kd + Kf) / Kp1)
    phi2P_ma = JP680_ma / (Q * a2_m)
    q2_ma = 1 - CB6F_ma / CB6F_m

    Kn2_ma = (
        (
            Kp2**2 * phi2P_ma**2
            - 2 * Kp2**2 * phi2P_ma * q2_ma
            + Kp2**2 * q2_ma**2
            - 4 * Kp2 * Ku2 * phi2P_ma**2 * q2_ma
            + 2 * Kp2 * Ku2 * phi2P_ma**2
            + 2 * Kp2 * Ku2 * phi2P_ma * q2_ma
            + Ku2**2 * phi2P_ma**2
        ) ** (1 / 2)
        - Kp2 * phi2P_ma
        + Ku2 * phi2P_ma
        + Kp2 * q2_ma
    ) / (2 * phi2P_ma) - Kf - Ku2 - Kd

    if _matlab_truth(Vpmax_m == 0):
        CB6F_s = 0
        RUB_s = 0
        CB6F_sa = 0
        phi1P_sa = 0
        q1_sa = 0
        phi2P_sa = 0
        q2_sa = 0
        Kn2_sa = 0

    if _matlab_truth(Vcmax_s > 0):
        CB6F_s = CB6F * vq_frac
        RUB_s = RUB * vc_frac
        CB6F_sa = JP700_sjj / kq
        phi1P_sa = JP700_sa / (Q * a1_s)
        q1_sa = phi1P_sa * ((Kp1 + Kd + Kf) / Kp1)
        phi2P_sa = JP680_sa / (Q * a2_s)
        q2_sa = 1 - CB6F_sa / CB6F_s

        Kn2_sa = (
            (
                Kp2**2 * phi2P_sa**2
                - 2 * Kp2**2 * phi2P_sa * q2_sa
                + Kp2**2 * q2_sa**2
                - 4 * Kp2 * Ku2 * phi2P_sa**2 * q2_sa
                + 2 * Kp2 * Ku2 * phi2P_sa**2
                + 2 * Kp2 * Ku2 * phi2P_sa * q2_sa
                + Ku2**2 * phi2P_sa**2
            ) ** (1 / 2)
            - Kp2 * phi2P_sa
            + Ku2 * phi2P_sa
            + Kp2 * q2_sa
        ) / (2 * phi2P_sa) - Kf - Ku2 - Kd

    phi2p_ma = q2_ma * Kp2 / (Kp2 + Kn2_ma + Kd + Kf + Ku2)
    phi2n_ma = q2_ma * Kn2_ma / (Kp2 + Kn2_ma + Kd + Kf + Ku2) + (1 - q2_ma) * Kn2_ma / (Kn2_ma + Kd + Kf + Ku2)
    phi2d_ma = q2_ma * Kd / (Kp2 + Kn2_ma + Kd + Kf + Ku2) + (1 - q2_ma) * Kd / (Kn2_ma + Kd + Kf + Ku2)
    phi2f_ma = q2_ma * Kf / (Kp2 + Kn2_ma + Kd + Kf + Ku2) + (1 - q2_ma) * Kf / (Kn2_ma + Kd + Kf + Ku2)
    phi2u_ma = q2_ma * Ku2 / (Kp2 + Kn2_ma + Kd + Kf + Ku2) + (1 - q2_ma) * Ku2 / (Kn2_ma + Kd + Kf + Ku2)
    phi2P_ma = phi2p_ma / (1 - phi2u_ma)
    phi2N_ma = phi2n_ma / (1 - phi2u_ma)
    phi2D_ma = phi2d_ma / (1 - phi2u_ma)
    phi2F_ma = phi2f_ma / (1 - phi2u_ma)

    phi1P_ma = q1_ma * Kp1 / (Kp1 + Kd + Kf)
    phi1N_ma = (1 - q1_ma) * Kn1 / (Kn1 + Kd + Kf)
    phi1D_ma = q1_ma * Kd / (Kp1 + Kd + Kf) + (1 - q1_ma) * Kd / (Kn1 + Kd + Kf)
    phi1F_ma = q1_ma * Kf / (Kp1 + Kd + Kf) + (1 - q1_ma) * Kf / (Kn1 + Kd + Kf)

    phi2p_sa = q2_sa * Kp2 / (Kp2 + Kn2_sa + Kd + Kf + Ku2)
    phi2n_sa = q2_sa * Kn2_sa / (Kp2 + Kn2_sa + Kd + Kf + Ku2) + (1 - q2_sa) * Kn2_sa / (Kn2_sa + Kd + Kf + Ku2)
    phi2d_sa = q2_sa * Kd / (Kp2 + Kn2_sa + Kd + Kf + Ku2) + (1 - q2_sa) * Kd / (Kn2_sa + Kd + Kf + Ku2)
    phi2f_sa = q2_sa * Kf / (Kp2 + Kn2_sa + Kd + Kf + Ku2) + (1 - q2_sa) * Kf / (Kn2_sa + Kd + Kf + Ku2)
    phi2u_sa = q2_sa * Ku2 / (Kp2 + Kn2_sa + Kd + Kf + Ku2) + (1 - q2_sa) * Ku2 / (Kn2_sa + Kd + Kf + Ku2)
    phi2P_sa = phi2p_sa / (1 - phi2u_sa)
    phi2N_sa = phi2n_sa / (1 - phi2u_sa)
    phi2D_sa = phi2d_sa / (1 - phi2u_sa)
    phi2F_sa = phi2f_sa / (1 - phi2u_sa)

    phi1P_sa = q1_sa * Kp1 / (Kp1 + Kd + Kf)
    phi1N_sa = (1 - q1_sa) * Kn1 / (Kn1 + Kd + Kf)
    phi1D_sa = q1_sa * Kd / (Kp1 + Kd + Kf) + (1 - q1_sa) * Kd / (Kn1 + Kd + Kf)
    phi1F_sa = q1_sa * Kf / (Kp1 + Kd + Kf) + (1 - q1_sa) * Kf / (Kn1 + Kd + Kf)

    Fs_ma = a2_m * phi2F_ma * eps2 + a1_m * phi1F_ma * eps1
    Fm_ma = a2_m * Kf / (Kd + Kf) * eps2 + a1_m * Kf / (Kn1 + Kd + Kf) * eps1
    Fo_ma = a2_m * Kf / (Kp2 + Kd + Kf) * eps2 + a1_m * Kf / (Kp1 + Kd + Kf) * eps1
    Fmp_ma = a2_m * Kf / (Kn2_ma + Kd + Kf) * eps2 + a1_m * Kf / (Kn1 + Kd + Kf) * eps1
    Fop_ma = a2_m * Kf / (Kp2 + Kn2_ma + Kd + Kf) * eps2 + a1_m * Kf / (Kp1 + Kd + Kf) * eps1

    Fs_sa = a2_s * phi2F_sa * eps2 + a1_s * phi1F_sa * eps1
    Fm_sa = a2_s * Kf / (Kd + Kf) * eps2 + a1_s * Kf / (Kn1 + Kd + Kf) * eps1
    Fo_sa = a2_s * Kf / (Kp2 + Kd + Kf) * eps2 + a1_s * Kf / (Kp1 + Kd + Kf) * eps1
    Fmp_sa = a2_s * Kf / (Kn2_sa + Kd + Kf) * eps2 + a1_s * Kf / (Kn1 + Kd + Kf) * eps1
    Fop_sa = a2_s * Kf / (Kp2 + Kn2_sa + Kd + Kf) * eps2 + a1_s * Kf / (Kp1 + Kd + Kf) * eps1

    Fs_a = Fs_ma + Fs_sa
    Fm_a = Fm_ma + Fm_sa
    Fo_a = Fo_ma + Fo_sa
    Fmp_a = Fmp_ma + Fmp_sa
    Fop_a = Fop_ma + Fop_sa

    PAM1_a = (Fmp_a - Fs_a) / (Fmp_a - Fop_a)
    PAM2_a = (Fmp_a - Fs_a) * Fop_a / ((Fmp_a - Fop_a) * Fs_a)
    PAM3_a = Fm_a / Fmp_a - 1
    PAM4_a = 1 - Fs_a / Fmp_a
    PAM5_a = Fs_a * (1 / Fmp_a - 1 / Fm_a)
    PAM6_a = Fs_a / Fm_a

    return workspace2struct_fun(
        exclude={
            "v",
            "Abs",
            "CB6F",
            "RUB",
            "Rdsc",
            "Kf",
            "Kd",
            "Kp1",
            "Kn1",
            "Kp2",
            "Ku2",
            "kq",
            "nl",
            "nc",
            "kc",
            "ko",
            "Kc",
            "Ko",
            "c3c4_solve_cc",
            "c3c4_solve_cj",
            "c3c4_solve_jc",
            "c3c4_solve_jj",
            "c4_solve_cc",
            "c4_solve_cj",
            "c4_solve_jc",
            "c4_solve_jj",
            "eps1",
            "eps2",
        }
    )
