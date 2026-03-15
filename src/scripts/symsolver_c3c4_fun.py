from datetime import datetime, timezone
from functools import lru_cache
from pathlib import Path
from types import SimpleNamespace

import cloudpickle
import numpy as np
import sympy as sp

_SYMSOLVER_CACHE_DIR = Path(__file__).resolve().parents[2] / "output" / "symsolver_cache"
_SYMSOLVER_CACHE_VERSION = 3


def _strcmp(pathway_option, target):
    if isinstance(pathway_option, (list, tuple)) and len(pathway_option) == 1:
        pathway_option = pathway_option[0]
    return pathway_option == target


def _normalize_pathway_option(pathway_option):
    if isinstance(pathway_option, (list, tuple)) and len(pathway_option) == 1:
        return pathway_option[0]
    return pathway_option


def _solve_first(equation, symbol):
    solutions = sp.solve(equation, symbol)
    if isinstance(solutions, dict):
        return solutions[symbol]
    return solutions[0]


def _lambdify(expression, args, real_output=False):
    function = sp.lambdify(args, expression, modules="numpy")

    if not real_output:
        return function

    def wrapped(*values):
        return np.real(function(*values))

    return wrapped


def _lambdify_complex_real(expression, args):
    function = sp.lambdify(args, expression, modules="numpy")

    def wrapped(*values):
        complex_values = [np.asarray(value, dtype=np.complex128) for value in values]
        return np.real(np.asarray(function(*complex_values), dtype=np.complex128))

    return wrapped


def _balance_polynomial(balance_expression, symbol):
    numerator, _ = sp.fraction(sp.cancel(sp.together(balance_expression)))
    return sp.Poly(sp.expand(numerator), symbol, domain="EX")


def _solve_balance_polynomial(balance_expression, symbol):
    try:
        polynomial = _balance_polynomial(balance_expression, symbol)
        return sp.solve(
            polynomial.as_expr(),
            symbol,
            simplify=False,
            check=False,
            cubics=True,
            quartics=True,
        )
    except sp.PolynomialError:
        numerator, _ = sp.fraction(sp.cancel(sp.together(balance_expression)))
        return sp.solve(
            sp.expand(numerator),
            symbol,
            simplify=False,
            check=False,
            cubics=True,
            quartics=True,
        )


def symsolver_c3c4_fun(pathway_option):
    return _load_or_solve_symsolver_c3c4_fun(_normalize_pathway_option(pathway_option))


@lru_cache(maxsize=None)
def _load_or_solve_symsolver_c3c4_fun(pathway_option):
    cached_namespace = load_symsolver_c3c4_cache(pathway_option)
    if cached_namespace is not None:
        return cached_namespace

    return _symsolver_c3c4_fun(pathway_option)


def _slugify_pathway_option(pathway_option):
    return (
        _normalize_pathway_option(pathway_option)
        .lower()
        .replace("-", "_")
        .replace(" ", "_")
    )


def get_symsolver_c3c4_cache_path(pathway_option):
    return _SYMSOLVER_CACHE_DIR / f"{_slugify_pathway_option(pathway_option)}.pkl"


def load_symsolver_c3c4_cache(pathway_option):
    cache_path = get_symsolver_c3c4_cache_path(pathway_option)
    if not cache_path.exists():
        return None

    with cache_path.open("rb") as handle:
        payload = cloudpickle.load(handle)

    if isinstance(payload, SimpleNamespace):
        return None

    if isinstance(payload, dict) and "namespace" in payload:
        if payload.get("cache_version") != _SYMSOLVER_CACHE_VERSION:
            return None
        return payload["namespace"]

    return None


def build_symsolver_c3c4_cache(pathway_option, overwrite=False):
    pathway_option = _normalize_pathway_option(pathway_option)
    cache_path = get_symsolver_c3c4_cache_path(pathway_option)

    if cache_path.exists() and not overwrite:
        try:
            with cache_path.open("rb") as handle:
                payload = cloudpickle.load(handle)
        except Exception:
            payload = None

        if (
            isinstance(payload, dict)
            and payload.get("cache_version") == _SYMSOLVER_CACHE_VERSION
        ):
            raise FileExistsError(f"Cache already exists: {cache_path}")

    namespace = _symsolver_c3c4_fun(pathway_option)

    cache_path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "cache_version": _SYMSOLVER_CACHE_VERSION,
        "pathway_option": pathway_option,
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "namespace": namespace,
    }
    with cache_path.open("wb") as handle:
        cloudpickle.dump(payload, handle)

    _load_or_solve_symsolver_c3c4_fun.cache_clear()
    return cache_path


@lru_cache(maxsize=None)
def _symsolver_c3c4_fun(pathway_option):
    (
        A_gross_s_c,
        A_net_s_c,
        A_gross_s_j,
        A_net_s_j,
        C_s,
        CO2_m,
        g_bs_CO2,
        g_bs_O2,
        J_PSI_s_j,
        J_PSI_m_j,
        J_PSII_s_j,
        J_PSII_m_j,
        K_m_CO2,
        K_m_O2,
        K_m_PEPC_CO2,
        k_P1,
        k_P2,
        k_D,
        k_F,
        k_U2,
        k_cat_CO2,
        k_cat_O2,
        k_q,
        ATP_e_ratio_linear,
        ATP_e_ratio_cyclic,
        O_s,
        O2_m,
        Pressure,
        phi_P2_m_actual,
        phi_p2_m_actual,
        phi_u2_m_actual,
        phi_P2_s_actual,
        phi_p2_s_actual,
        phi_u2_s_actual,
        R_d_s,
        Specificity,
        Vg_c,
        Vp_c,
        Vg_j,
        Vp_j,
        V_c_m_c,
        V_o_m_c,
        V_c_s_c,
        V_o_s_c,
        V_c_m_j,
        V_o_m_j,
        V_c_s_j,
        V_o_s_j,
        V_p_max_m,
        V_c_max_s,
        V_c_max_m,
    ) = sp.symbols(
        "A_gross_s_c A_net_s_c A_gross_s_j A_net_s_j C_s CO2_m g_bs_CO2 g_bs_O2 J_PSI_s_j J_PSI_m_j "
        "J_PSII_s_j J_PSII_m_j K_m_CO2 K_m_O2 K_m_PEPC_CO2 k_P1 k_P2 k_D k_F k_U2 k_cat_CO2 k_cat_O2 k_q ATP_e_ratio_linear ATP_e_ratio_cyclic "
        "O_s O2_m Pressure phi_P2_m_actual phi_p2_m_actual phi_u2_m_actual phi_P2_s_actual phi_p2_s_actual phi_u2_s_actual "
        "R_d_s Specificity Vg_c Vp_c Vg_j Vp_j V_c_m_c V_o_m_c V_c_s_c V_o_s_c V_c_m_j V_o_m_j "
        "V_c_s_j V_o_s_j V_p_max_m V_c_max_s V_c_max_m"
    )

    if _strcmp(pathway_option, "Type-I-C3-C4"):
        Vg_j_expr = J_PSII_m_j * O2_m / (2 * Specificity * CO2_m) / (4 * (1 + O2_m / (Specificity * CO2_m)))
        Vg_c_expr = V_c_max_m * O2_m / (2 * Specificity) / (CO2_m + K_m_CO2 * (1 + O2_m / K_m_O2))

        CO2_s_jj = CO2_m + (Vg_j_expr - A_net_s_j) * Pressure / g_bs_CO2
        O2_s_jj = O2_m + A_net_s_j * Pressure / g_bs_O2
        ratio_sjj = O2_s_jj / (Specificity * CO2_s_jj)
        J_PSII_s_jj = J_PSI_s_j / (
            1
            - (ATP_e_ratio_linear / ATP_e_ratio_cyclic)
            + (3 + sp.Rational(7, 2) * ratio_sjj) / ((4 + 4 * ratio_sjj) * ATP_e_ratio_cyclic)
        )
        V_c_s_jj = J_PSII_s_jj / (4 * (1 + ratio_sjj))
        A_gross_s_jj = V_c_s_jj * (1 - ratio_sjj / 2)
        c3c4_soln_jj = _solve_balance_polynomial(A_gross_s_jj - R_d_s - A_net_s_j, A_net_s_j)

        CO2_s_jc = CO2_m + (Vg_j_expr - A_net_s_c) * Pressure / g_bs_CO2
        O2_s_jc = O2_m + A_net_s_c * Pressure / g_bs_O2
        V_c_s_jc = CO2_s_jc * V_c_max_s / (CO2_s_jc + K_m_CO2 * (1 + O2_s_jc / K_m_O2))
        A_gross_s_jc = V_c_s_jc * (1 - O2_s_jc / (2 * Specificity * CO2_s_jc))
        c3c4_soln_jc = _solve_balance_polynomial(A_gross_s_jc - R_d_s - A_net_s_c, A_net_s_c)

        CO2_s_cj = CO2_m + (Vg_c_expr - A_net_s_j) * Pressure / g_bs_CO2
        O2_s_cj = O2_m + A_net_s_j * Pressure / g_bs_O2
        ratio_scj = O2_s_cj / (Specificity * CO2_s_cj)
        J_PSII_s_cj = J_PSI_s_j / (
            1
            - (ATP_e_ratio_linear / ATP_e_ratio_cyclic)
            + (3 + sp.Rational(7, 2) * ratio_scj) / ((4 + 4 * ratio_scj) * ATP_e_ratio_cyclic)
        )
        V_c_s_cj = J_PSII_s_cj / (4 * (1 + ratio_scj))
        A_gross_s_cj = V_c_s_cj * (1 - ratio_scj / 2)
        c3c4_soln_cj = _solve_balance_polynomial(A_gross_s_cj - R_d_s - A_net_s_j, A_net_s_j)

        CO2_s_cc = CO2_m + (Vg_c_expr - A_net_s_c) * Pressure / g_bs_CO2
        O2_s_cc = O2_m + A_net_s_c * Pressure / g_bs_O2
        V_c_s_cc = CO2_s_cc * V_c_max_s / (CO2_s_cc + K_m_CO2 * (1 + O2_s_cc / K_m_O2))
        A_gross_s_cc = V_c_s_cc * (1 - O2_s_cc / (2 * Specificity * CO2_s_cc))
        c3c4_soln_cc = _solve_balance_polynomial(A_gross_s_cc - R_d_s - A_net_s_c, A_net_s_c)

        solve_C3C4_jj = _lambdify(
            c3c4_soln_jj[1],
            (CO2_m, J_PSII_m_j, J_PSI_s_j, O2_m, Pressure, R_d_s, Specificity, g_bs_CO2, g_bs_O2, ATP_e_ratio_cyclic, ATP_e_ratio_linear),
        )
        solve_C3C4_jc = _lambdify(
            c3c4_soln_jc[0],
            (CO2_m, J_PSII_m_j, K_m_CO2, K_m_O2, O2_m, Pressure, R_d_s, Specificity, V_c_max_s, g_bs_CO2, g_bs_O2),
        )
        solve_C3C4_cj = _lambdify(
            c3c4_soln_cj[1],
            (CO2_m, J_PSI_s_j, K_m_CO2, K_m_O2, O2_m, Pressure, R_d_s, Specificity, V_c_max_m, g_bs_CO2, g_bs_O2, ATP_e_ratio_cyclic, ATP_e_ratio_linear),
        )
        solve_C3C4_cc = _lambdify(
            c3c4_soln_cc[0],
            (CO2_m, K_m_CO2, K_m_O2, O2_m, Pressure, R_d_s, Specificity, V_c_max_m, V_c_max_s, g_bs_CO2, g_bs_O2),
        )

        return SimpleNamespace(
            solve_C3C4_jj=solve_C3C4_jj,
            solve_C3C4_jc=solve_C3C4_jc,
            solve_C3C4_cj=solve_C3C4_cj,
            solve_C3C4_cc=solve_C3C4_cc,
        )

    if _strcmp(pathway_option, "NADP-ME-C4"):
        Vp_j_expr = J_PSI_m_j / (2 * (1 - (ATP_e_ratio_linear / ATP_e_ratio_cyclic) + 1 / ATP_e_ratio_cyclic))
        Vp_c_expr = V_p_max_m * CO2_m / (K_m_PEPC_CO2 + CO2_m)

        CO2_s_jj = CO2_m + (Vp_j_expr - A_net_s_j) * Pressure / g_bs_CO2
        O2_s_jj = O2_m + (A_net_s_j - Vp_j_expr / 2) * Pressure / g_bs_O2
        ratio_sjj = O2_s_jj / (Specificity * CO2_s_jj)
        alpha_sjj = (3 + sp.Rational(7, 2) * ratio_sjj) / ((4 + 4 * ratio_sjj) * ATP_e_ratio_cyclic)
        J_PSII_s_jj = (J_PSI_s_j - (Vp_j_expr / 2) * alpha_sjj) / (
            1 - (ATP_e_ratio_linear / ATP_e_ratio_cyclic) + alpha_sjj
        )
        V_c_s_jj = (J_PSII_s_jj + 2 * Vp_j_expr) / (4 * (1 + ratio_sjj))
        A_gross_s_jj = V_c_s_jj * (1 - ratio_sjj / 2)
        c4_soln_jj = _solve_balance_polynomial(A_gross_s_jj - R_d_s - A_net_s_j, A_net_s_j)

        CO2_s_jc = CO2_m + (Vp_j_expr - A_net_s_c) * Pressure / g_bs_CO2
        O2_s_jc = O2_m + (A_net_s_c - Vp_j_expr / 2) * Pressure / g_bs_O2
        V_c_s_jc = CO2_s_jc * V_c_max_s / (CO2_s_jc + K_m_CO2 * (1 + O2_s_jc / K_m_O2))
        A_gross_s_jc = V_c_s_jc * (1 - O2_s_jc / (2 * Specificity * CO2_s_jc))
        c4_soln_jc = _solve_balance_polynomial(A_gross_s_jc - R_d_s - A_net_s_c, A_net_s_c)

        CO2_s_cj = CO2_m + (Vp_c_expr - A_net_s_j) * Pressure / g_bs_CO2
        O2_s_cj = O2_m + (A_net_s_j - Vp_c_expr / 2) * Pressure / g_bs_O2
        ratio_scj = O2_s_cj / (Specificity * CO2_s_cj)
        alpha_scj = (3 + sp.Rational(7, 2) * ratio_scj) / ((4 + 4 * ratio_scj) * ATP_e_ratio_cyclic)
        J_PSII_s_cj = (J_PSI_s_j - (Vp_c_expr / 2) * alpha_scj) / (
            1 - (ATP_e_ratio_linear / ATP_e_ratio_cyclic) + alpha_scj
        )
        V_c_s_cj = (J_PSII_s_cj + 2 * Vp_c_expr) / (4 * (1 + ratio_scj))
        A_gross_s_cj = V_c_s_cj * (1 - ratio_scj / 2)
        c4_soln_cj = _solve_balance_polynomial(A_gross_s_cj - R_d_s - A_net_s_j, A_net_s_j)

        CO2_s_cc = CO2_m + (Vp_c_expr - A_net_s_c) * Pressure / g_bs_CO2
        O2_s_cc = O2_m + (A_net_s_c - Vp_c_expr / 2) * Pressure / g_bs_O2
        V_c_s_cc = CO2_s_cc * V_c_max_s / (CO2_s_cc + K_m_CO2 * (1 + O2_s_cc / K_m_O2))
        A_gross_s_cc = V_c_s_cc * (1 - O2_s_cc / (2 * Specificity * CO2_s_cc))
        c4_soln_cc = _solve_balance_polynomial(A_gross_s_cc - R_d_s - A_net_s_c, A_net_s_c)

        solve_C4_jj = _lambdify_complex_real(
            c4_soln_jj[0],
            (CO2_m, J_PSI_m_j, J_PSI_s_j, O2_m, Pressure, R_d_s, Specificity, g_bs_CO2, g_bs_O2, ATP_e_ratio_cyclic, ATP_e_ratio_linear),
        )
        solve_C4_jc = _lambdify(
            c4_soln_jc[0],
            (CO2_m, J_PSI_m_j, K_m_CO2, K_m_O2, O2_m, Pressure, R_d_s, Specificity, V_c_max_s, g_bs_CO2, g_bs_O2, ATP_e_ratio_cyclic, ATP_e_ratio_linear),
        )
        solve_C4_cj = _lambdify_complex_real(
            c4_soln_cj[0],
            (CO2_m, J_PSI_s_j, K_m_PEPC_CO2, O2_m, Pressure, R_d_s, Specificity, V_p_max_m, g_bs_CO2, g_bs_O2, ATP_e_ratio_cyclic, ATP_e_ratio_linear),
        )
        solve_C4_cc = _lambdify(
            c4_soln_cc[0],
            (CO2_m, K_m_CO2, K_m_O2, K_m_PEPC_CO2, O2_m, Pressure, R_d_s, Specificity, V_c_max_s, V_p_max_m, g_bs_CO2, g_bs_O2),
        )

        return SimpleNamespace(
            solve_C4_jj=solve_C4_jj,
            solve_C4_jc=solve_C4_jc,
            solve_C4_cj=solve_C4_cj,
            solve_C4_cc=solve_C4_cc,
        )

    return SimpleNamespace()
