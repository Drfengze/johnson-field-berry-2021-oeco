from datetime import datetime, timezone
from functools import lru_cache
from pathlib import Path
from types import SimpleNamespace

import cloudpickle
import numpy as np
import sympy as sp

_SYMSOLVER_CACHE_DIR = Path(__file__).resolve().parents[2] / "output" / "symsolver_cache"
_SYMSOLVER_CACHE_VERSION = 2


def _strcmp(pathway_opt, target):
    if isinstance(pathway_opt, (list, tuple)) and len(pathway_opt) == 1:
        pathway_opt = pathway_opt[0]
    return pathway_opt == target


def _normalize_pathway_opt(pathway_opt):
    if isinstance(pathway_opt, (list, tuple)) and len(pathway_opt) == 1:
        return pathway_opt[0]
    return pathway_opt


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


def symsolver_c3c4_fun(pathway_opt):
    return _load_or_solve_symsolver_c3c4_fun(_normalize_pathway_opt(pathway_opt))


@lru_cache(maxsize=None)
def _load_or_solve_symsolver_c3c4_fun(pathway_opt):
    cached_namespace = load_symsolver_c3c4_cache(pathway_opt)
    if cached_namespace is not None:
        return cached_namespace

    return _symsolver_c3c4_fun(pathway_opt)


def _slugify_pathway_opt(pathway_opt):
    return (
        _normalize_pathway_opt(pathway_opt)
        .lower()
        .replace("-", "_")
        .replace(" ", "_")
    )


def get_symsolver_c3c4_cache_path(pathway_opt):
    return _SYMSOLVER_CACHE_DIR / f"{_slugify_pathway_opt(pathway_opt)}.pkl"


def load_symsolver_c3c4_cache(pathway_opt):
    cache_path = get_symsolver_c3c4_cache_path(pathway_opt)
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


def build_symsolver_c3c4_cache(pathway_opt, overwrite=False):
    pathway_opt = _normalize_pathway_opt(pathway_opt)
    cache_path = get_symsolver_c3c4_cache_path(pathway_opt)

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

    namespace = _symsolver_c3c4_fun(pathway_opt)

    cache_path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "cache_version": _SYMSOLVER_CACHE_VERSION,
        "pathway_opt": pathway_opt,
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "namespace": namespace,
    }
    with cache_path.open("wb") as handle:
        cloudpickle.dump(payload, handle)

    _load_or_solve_symsolver_c3c4_fun.cache_clear()
    return cache_path


@lru_cache(maxsize=None)
def _symsolver_c3c4_fun(pathway_opt):
    (
        Ag_sc,
        An_sc,
        Ag_sj,
        An_sj,
        C_s,
        C_m,
        gbs,
        gbso,
        JP700_sj,
        JP700_mj,
        JP680_sj,
        JP680_mj,
        Kc,
        Ko,
        Kp,
        Kp1,
        Kp2,
        Kd,
        Kf,
        Ku2,
        kc,
        ko,
        kq,
        nl,
        nc,
        O_s,
        O_m,
        P,
        phi2P_ma,
        phi2p_ma,
        phi2u_ma,
        phi2P_sa,
        phi2p_sa,
        phi2u_sa,
        Rd_s,
        S,
        Vg_c,
        Vp_c,
        Vg_j,
        Vp_j,
        Vc_mc,
        Vo_mc,
        Vc_sc,
        Vo_sc,
        Vc_mj,
        Vo_mj,
        Vc_sj,
        Vo_sj,
        Vpmax_m,
        Vcmax_s,
        Vcmax_m,
    ) = sp.symbols(
        "Ag_sc An_sc Ag_sj An_sj C_s C_m gbs gbso JP700_sj JP700_mj "
        "JP680_sj JP680_mj Kc Ko Kp Kp1 Kp2 Kd Kf Ku2 kc ko kq nl nc "
        "O_s O_m P phi2P_ma phi2p_ma phi2u_ma phi2P_sa phi2p_sa phi2u_sa "
        "Rd_s S Vg_c Vp_c Vg_j Vp_j Vc_mc Vo_mc Vc_sc Vo_sc Vc_mj Vo_mj "
        "Vc_sj Vo_sj Vpmax_m Vcmax_s Vcmax_m"
    )

    if _strcmp(pathway_opt, "Type-I-C3-C4"):
        Vg_j_expr = JP680_mj * O_m / (2 * S * C_m) / (4 * (1 + O_m / (S * C_m)))
        Vg_c_expr = Vcmax_m * O_m / (2 * S) / (C_m + Kc * (1 + O_m / Ko))

        C_sjj = C_m + (Vg_j_expr - An_sj) * P / gbs
        O_sjj = O_m + An_sj * P / gbso
        ratio_sjj = O_sjj / (S * C_sjj)
        JP680_sjj = JP700_sj / (
            1
            - (nl / nc)
            + (3 + sp.Rational(7, 2) * ratio_sjj) / ((4 + 4 * ratio_sjj) * nc)
        )
        Vc_sjj = JP680_sjj / (4 * (1 + ratio_sjj))
        Ag_sjj = Vc_sjj * (1 - ratio_sjj / 2)
        c3c4_soln_jj = _solve_balance_polynomial(Ag_sjj - Rd_s - An_sj, An_sj)

        C_sjc = C_m + (Vg_j_expr - An_sc) * P / gbs
        O_sjc = O_m + An_sc * P / gbso
        Vc_sjc = C_sjc * Vcmax_s / (C_sjc + Kc * (1 + O_sjc / Ko))
        Ag_sjc = Vc_sjc * (1 - O_sjc / (2 * S * C_sjc))
        c3c4_soln_jc = _solve_balance_polynomial(Ag_sjc - Rd_s - An_sc, An_sc)

        C_scj = C_m + (Vg_c_expr - An_sj) * P / gbs
        O_scj = O_m + An_sj * P / gbso
        ratio_scj = O_scj / (S * C_scj)
        JP680_scj = JP700_sj / (
            1
            - (nl / nc)
            + (3 + sp.Rational(7, 2) * ratio_scj) / ((4 + 4 * ratio_scj) * nc)
        )
        Vc_scj = JP680_scj / (4 * (1 + ratio_scj))
        Ag_scj = Vc_scj * (1 - ratio_scj / 2)
        c3c4_soln_cj = _solve_balance_polynomial(Ag_scj - Rd_s - An_sj, An_sj)

        C_scc = C_m + (Vg_c_expr - An_sc) * P / gbs
        O_scc = O_m + An_sc * P / gbso
        Vc_scc = C_scc * Vcmax_s / (C_scc + Kc * (1 + O_scc / Ko))
        Ag_scc = Vc_scc * (1 - O_scc / (2 * S * C_scc))
        c3c4_soln_cc = _solve_balance_polynomial(Ag_scc - Rd_s - An_sc, An_sc)

        c3c4_solve_jj = _lambdify(
            c3c4_soln_jj[1],
            (C_m, JP680_mj, JP700_sj, O_m, P, Rd_s, S, gbs, gbso, nc, nl),
        )
        c3c4_solve_jc = _lambdify(
            c3c4_soln_jc[0],
            (C_m, JP680_mj, Kc, Ko, O_m, P, Rd_s, S, Vcmax_s, gbs, gbso),
        )
        c3c4_solve_cj = _lambdify(
            c3c4_soln_cj[1],
            (C_m, JP700_sj, Kc, Ko, O_m, P, Rd_s, S, Vcmax_m, gbs, gbso, nc, nl),
        )
        c3c4_solve_cc = _lambdify(
            c3c4_soln_cc[0],
            (C_m, Kc, Ko, O_m, P, Rd_s, S, Vcmax_m, Vcmax_s, gbs, gbso),
        )

        return SimpleNamespace(
            c3c4_solve_jj=c3c4_solve_jj,
            c3c4_solve_jc=c3c4_solve_jc,
            c3c4_solve_cj=c3c4_solve_cj,
            c3c4_solve_cc=c3c4_solve_cc,
        )

    if _strcmp(pathway_opt, "NADP-ME-C4"):
        Vp_j_expr = JP700_mj / (2 * (1 - (nl / nc) + 1 / nc))
        Vp_c_expr = Vpmax_m * C_m / (Kp + C_m)

        C_sjj = C_m + (Vp_j_expr - An_sj) * P / gbs
        O_sjj = O_m + (An_sj - Vp_j_expr / 2) * P / gbso
        ratio_sjj = O_sjj / (S * C_sjj)
        alpha_sjj = (3 + sp.Rational(7, 2) * ratio_sjj) / ((4 + 4 * ratio_sjj) * nc)
        JP680_sjj = (JP700_sj - (Vp_j_expr / 2) * alpha_sjj) / (
            1 - (nl / nc) + alpha_sjj
        )
        Vc_sjj = (JP680_sjj + 2 * Vp_j_expr) / (4 * (1 + ratio_sjj))
        Ag_sjj = Vc_sjj * (1 - ratio_sjj / 2)
        c4_soln_jj = _solve_balance_polynomial(Ag_sjj - Rd_s - An_sj, An_sj)

        C_sjc = C_m + (Vp_j_expr - An_sc) * P / gbs
        O_sjc = O_m + (An_sc - Vp_j_expr / 2) * P / gbso
        Vc_sjc = C_sjc * Vcmax_s / (C_sjc + Kc * (1 + O_sjc / Ko))
        Ag_sjc = Vc_sjc * (1 - O_sjc / (2 * S * C_sjc))
        c4_soln_jc = _solve_balance_polynomial(Ag_sjc - Rd_s - An_sc, An_sc)

        C_scj = C_m + (Vp_c_expr - An_sj) * P / gbs
        O_scj = O_m + (An_sj - Vp_c_expr / 2) * P / gbso
        ratio_scj = O_scj / (S * C_scj)
        alpha_scj = (3 + sp.Rational(7, 2) * ratio_scj) / ((4 + 4 * ratio_scj) * nc)
        JP680_scj = (JP700_sj - (Vp_c_expr / 2) * alpha_scj) / (
            1 - (nl / nc) + alpha_scj
        )
        Vc_scj = (JP680_scj + 2 * Vp_c_expr) / (4 * (1 + ratio_scj))
        Ag_scj = Vc_scj * (1 - ratio_scj / 2)
        c4_soln_cj = _solve_balance_polynomial(Ag_scj - Rd_s - An_sj, An_sj)

        C_scc = C_m + (Vp_c_expr - An_sc) * P / gbs
        O_scc = O_m + (An_sc - Vp_c_expr / 2) * P / gbso
        Vc_scc = C_scc * Vcmax_s / (C_scc + Kc * (1 + O_scc / Ko))
        Ag_scc = Vc_scc * (1 - O_scc / (2 * S * C_scc))
        c4_soln_cc = _solve_balance_polynomial(Ag_scc - Rd_s - An_sc, An_sc)

        c4_solve_jj = _lambdify_complex_real(
            c4_soln_jj[0],
            (C_m, JP700_mj, JP700_sj, O_m, P, Rd_s, S, gbs, gbso, nc, nl),
        )
        c4_solve_jc = _lambdify(
            c4_soln_jc[0],
            (C_m, JP700_mj, Kc, Ko, O_m, P, Rd_s, S, Vcmax_s, gbs, gbso, nc, nl),
        )
        c4_solve_cj = _lambdify_complex_real(
            c4_soln_cj[0],
            (C_m, JP700_sj, Kp, O_m, P, Rd_s, S, Vpmax_m, gbs, gbso, nc, nl),
        )
        c4_solve_cc = _lambdify(
            c4_soln_cc[0],
            (C_m, Kc, Ko, Kp, O_m, P, Rd_s, S, Vcmax_s, Vpmax_m, gbs, gbso),
        )

        return SimpleNamespace(
            c4_solve_jj=c4_solve_jj,
            c4_solve_jc=c4_solve_jc,
            c4_solve_cj=c4_solve_cj,
            c4_solve_cc=c4_solve_cc,
        )

    return SimpleNamespace()
