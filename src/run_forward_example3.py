#!/usr/bin/env python3

import argparse
import time
from pathlib import Path
from types import SimpleNamespace

import matplotlib
import numpy as np
from scipy.io import savemat

matplotlib.use("Agg")

from scripts.configure_fun import configure_fun
from scripts.model_fun_c3c4 import model_fun_c3c4
from scripts.plotter_forward_fun_c3c4 import plotter_forward_fun_c3c4
from scripts.symsolver_c3c4_fun import symsolver_c3c4_fun

REPO_DIR = Path(__file__).resolve().parents[1]
OUTPUT_DIR = REPO_DIR / "output"
SUPPORTED_PATHWAYS = ("C3", "Type-I-C3-C4", "NADP-ME-C4")


def _parse_args():
    parser = argparse.ArgumentParser(
        description="Python translation of run_forward_example3.m."
    )
    parser.add_argument(
        "--pathway",
        choices=SUPPORTED_PATHWAYS,
        default="Type-I-C3-C4",
        help="Pathway option to simulate.",
    )
    parser.add_argument(
        "--output-suffix",
        default="",
        help="Optional suffix appended to the output directory name.",
    )
    return parser.parse_args()


def _build_example3_data(n):
    return SimpleNamespace(
        Tin=np.linspace(10, 40, n),
        Qin=np.full(n, 1200.0),
        Cin=np.full(n, 250.0),
        Oin=np.full(n, 0.209),
        Pin=np.full(n, 1.0),
    )


def _serialize_for_mat(value):
    if isinstance(value, SimpleNamespace):
        return {name: _serialize_for_mat(field_value) for name, field_value in vars(value).items()}

    if isinstance(value, dict):
        return {name: _serialize_for_mat(field_value) for name, field_value in value.items()}

    if callable(value):
        return getattr(value, "__name__", "python_function")

    if isinstance(value, (list, tuple)) and len(value) == 1 and isinstance(value[0], str):
        return value[0]

    return value

def main():
    args = _parse_args()
    start = time.time()

    currdir = REPO_DIR
    resultsdir = OUTPUT_DIR

    pathway_opt = args.pathway
    outputname = f"Example3-Temperature-{pathway_opt}"
    outputdir = resultsdir / f"{outputname}{args.output_suffix}"
    outputdir.mkdir(parents=True, exist_ok=True)

    n = 100
    data = _build_example3_data(n)

    v = configure_fun(data)
    v.Ku2 = 2e09

    if pathway_opt == "C3":
        v.CB6F = 175.0 / v.kq * 1e-06
        v.RUB = 50.0 / v.kc * 1e-06
        v.abs_frac = 0.0
        v.vq_frac = 0.0
        v.vc_frac = 0.0
        v.a2_s_frac = 0.0
        v.a2_m_frac = 0.52

    if pathway_opt == "Type-I-C3-C4":
        ss = symsolver_c3c4_fun(pathway_opt)
        v.CB6F = 175.0 / v.kq * 1e-06
        v.RUB = 50.0 / v.kc * 1e-06
        v.abs_frac = 0.05
        v.vq_frac = 0.05
        v.vc_frac = 0.1
        v.a2_s_frac = 0.52
        v.a2_m_frac = 0.52

    if pathway_opt == "NADP-ME-C4":
        ss = symsolver_c3c4_fun(pathway_opt)
        v.CB6F = 175.0 / v.kq * 1e-06
        v.RUB = 30.0 / v.kc * 1e-06
        v.Vpmax = v.RUB * v.kc * 2
        v.abs_frac = 0.4
        v.vq_frac = 0.4
        v.vc_frac = 1
        v.a2_s_frac = 0.47
        v.a2_m_frac = 0.47

    v.Model_id = "model_fun_c3c4"
    v.pathway_opt = pathway_opt

    if pathway_opt == "Type-I-C3-C4":
        v.c3c4_solve_cc = ss.c3c4_solve_cc
        v.c3c4_solve_cj = ss.c3c4_solve_cj
        v.c3c4_solve_jc = ss.c3c4_solve_jc
        v.c3c4_solve_jj = ss.c3c4_solve_jj

    if pathway_opt == "NADP-ME-C4":
        v.c4_solve_cc = ss.c4_solve_cc
        v.c4_solve_cj = ss.c4_solve_cj
        v.c4_solve_jc = ss.c4_solve_jc
        v.c4_solve_jj = ss.c4_solve_jj

    m = model_fun_c3c4(v)
    s = plotter_forward_fun_c3c4(outputname, v, m)

    s.figure1.set_size_inches(12, 12)
    s.figure1.savefig(outputdir / f"{outputname}-figure.png", dpi=300)

    savemat(outputdir / f"{outputname}-modelinputs.mat", {"v": _serialize_for_mat(v)})
    savemat(outputdir / f"{outputname}-modeloutputs.mat", {"m": _serialize_for_mat(m)})

    print(f"[done] output_dir={outputdir}")
    print(f"[done] elapsed_seconds={time.time() - start:.3f}")
    _ = currdir
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
