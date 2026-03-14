from types import SimpleNamespace

import matplotlib.pyplot as plt
import numpy as np


def _field(value, name):
    if isinstance(value, dict):
        return value[name]
    return getattr(value, name)


def _strcmp(pathway_opt, target):
    if isinstance(pathway_opt, (list, tuple)) and len(pathway_opt) == 1:
        pathway_opt = pathway_opt[0]
    return pathway_opt == target


def _array(value):
    return np.asarray(value, dtype=float)


def _series(x, value):
    value = np.asarray(value, dtype=float)
    if value.ndim == 0 or value.size == 1:
        return np.full(np.asarray(x).shape, float(value))
    return value


def _annotate(ax, xpos, ypos, label):
    xlim_curr = ax.get_xlim()
    ylim_curr = ax.get_ylim()
    ax.text(
        xlim_curr[0] + (xlim_curr[1] - xlim_curr[0]) * xpos,
        ylim_curr[0] + (ylim_curr[1] - ylim_curr[0]) * ypos,
        label,
    )


def plotter_forward_fun_c3c4(outputname, v, m):
    if np.ptp(_array(_field(m, "Q"))) > 0:
        x = _array(_field(m, "Q")) * 1e6
        xlab = "PAR (umol PPFD m-2 s-1)"

    if np.ptp(_array(_field(m, "C_m"))) > 0:
        x = _array(_field(m, "C_m")) * 1e6
        xlab = "Cm (ubar CO2)"

    if np.ptp(_array(_field(m, "T"))) > 0:
        x = _array(_field(m, "T"))
        xlab = "Tleaf (C)"

    if np.ptp(_array(_field(m, "Q"))) > 0 and np.ptp(_array(_field(m, "T"))) > 0:
        x = np.linspace(1, 1440, 1440) / 60
        xlab = "Time (hours)"

    figure1 = plt.figure()
    xpos = 0.85
    ypos = 0.9

    ax = figure1.add_subplot(5, 4, 1)
    ax.plot(x, _series(x, _array(_field(m, "JP680_ma")) * 1e6), "-b")
    ax.plot(x, _series(x, _array(_field(m, "JP680_sa")) * 1e6), "-r")
    ax.set_xlim(np.min(x), np.max(x))
    ax.set_ylim(0, 250)
    ax.set_ylabel("PS2 ETR (umol e- m-2 s-1)")
    ax.legend(["Mesophyll", "Bundle sheath"])
    _annotate(ax, xpos, ypos, "(a)")

    ax = figure1.add_subplot(5, 4, 2)
    ax.plot(x, _series(x, _array(_field(m, "phi2P_ma")) * 100), "-b")
    ax.plot(x, _series(x, _array(_field(m, "phi2P_sa")) * 100), "-r")
    ax.set_xlim(np.min(x), np.max(x))
    ax.set_ylim(0, 100)
    ax.set_ylabel("Phi2P (%)")
    _annotate(ax, xpos, ypos, "(b)")

    ax = figure1.add_subplot(5, 4, 3)
    ax.plot(x, _series(x, _array(_field(m, "phi2N_ma")) * 100), "-b")
    ax.plot(x, _series(x, _array(_field(m, "phi2N_sa")) * 100), "-r")
    ax.set_xlim(np.min(x), np.max(x))
    ax.set_ylim(0, 100)
    ax.set_ylabel("Phi2N (%)")
    _annotate(ax, xpos, ypos, "(c)")

    ax = figure1.add_subplot(5, 4, 4)
    ax.plot(x, _series(x, (_array(_field(m, "phi2D_ma")) + _array(_field(m, "phi2F_ma"))) * 100), "-b")
    ax.plot(x, _series(x, (_array(_field(m, "phi2D_sa")) + _array(_field(m, "phi2F_sa"))) * 100), "-r")
    ax.set_xlim(np.min(x), np.max(x))
    ax.set_ylim(0, 100)
    ax.set_ylabel("Phi2DF (%)")
    _annotate(ax, xpos, ypos, "(d)")

    ax = figure1.add_subplot(5, 4, 5)
    ax.plot(x, _series(x, _array(_field(m, "JP700_ma")) * 1e6), "-b")
    ax.plot(x, _series(x, _array(_field(m, "JP700_sa")) * 1e6), "-r")
    ax.set_xlim(np.min(x), np.max(x))
    ax.set_ylim(0, 250)
    ax.set_ylabel("PS1 ETR (umol e- m-2 s-1)")
    _annotate(ax, xpos, ypos, "(e)")

    ax = figure1.add_subplot(5, 4, 6)
    ax.plot(x, _series(x, _array(_field(m, "phi1P_ma")) * 100), "-b")
    ax.plot(x, _series(x, _array(_field(m, "phi1P_sa")) * 100), "-r")
    ax.set_xlim(np.min(x), np.max(x))
    ax.set_ylim(0, 100)
    ax.set_ylabel("Phi1P (%)")
    _annotate(ax, xpos, ypos, "(f)")

    ax = figure1.add_subplot(5, 4, 7)
    ax.plot(x, _series(x, _array(_field(m, "phi1N_ma")) * 100), "-b")
    ax.plot(x, _series(x, _array(_field(m, "phi1N_sa")) * 100), "-r")
    ax.set_xlim(np.min(x), np.max(x))
    ax.set_ylim(0, 100)
    ax.set_ylabel("Phi1N (%)")
    _annotate(ax, xpos, ypos, "(g)")

    ax = figure1.add_subplot(5, 4, 8)
    ax.plot(x, _series(x, (_array(_field(m, "phi1D_ma")) + _array(_field(m, "phi1F_ma"))) * 100), "-b")
    ax.plot(x, _series(x, (_array(_field(m, "phi1D_sa")) + _array(_field(m, "phi1F_sa"))) * 100), "-r")
    ax.set_xlim(np.min(x), np.max(x))
    ax.set_ylim(0, 100)
    ax.set_ylabel("Phi1DF (%)")
    _annotate(ax, xpos, ypos, "(h)")

    ax = figure1.add_subplot(5, 4, 9)
    ax.plot(x, _series(x, _array(_field(m, "JP700_ma")) * 1e6), "-b")
    ax.plot(x, _series(x, _array(_field(m, "JP700_sa")) * 1e6), "-r")
    ax.set_xlim(np.min(x), np.max(x))
    ax.set_ylim(0, 250)
    ax.set_ylabel("Cyt b6f ETR (umol e- m-2 s-1)")
    _annotate(ax, xpos, ypos, "(i)")

    ax = figure1.add_subplot(5, 4, 10)
    ax.plot(x, _series(x, (1 - _array(_field(m, "q2_ma"))) * 100), "-b")
    ax.plot(x, _series(x, (1 - _array(_field(m, "q2_sa"))) * 100), "-r")
    ax.set_xlim(np.min(x), np.max(x))
    ax.set_ylim(0, 100)
    ax.set_ylabel("PQH2/[PQ+PQH2] (%)")
    _annotate(ax, xpos, ypos, "(j)")

    ax = figure1.add_subplot(5, 4, 11)
    ax.plot(x, _series(x, _array(_field(m, "JP700_ma")) / _array(_field(m, "JP700_mj")) * (_array(_field(m, "Vqmax_m")) / _array(_field(m, "CB6F_m")))), "-b")
    ax.plot(
        x,
        _series(
            x,
            _array(_field(m, "JP700_sa"))
            / (
                _array(_field(m, "JP700_sjj")) * (_array(_field(m, "which_JP700_ma")) == 1)
                + _array(_field(m, "JP700_scj")) * (_array(_field(m, "which_JP700_ma")) == 2)
            )
            * (_array(_field(m, "Vqmax_s")) / _array(_field(m, "CB6F_s"))),
        ),
        "-r",
    )
    ax.set_xlim(np.min(x), np.max(x))
    ax.set_ylim(0, 600)
    ax.set_ylabel("Cyt b6f rate constant (s-1)")
    _annotate(ax, xpos, ypos, "(k)")

    ax = figure1.add_subplot(5, 4, 12)
    ax.plot(x, _series(x, (1 - _array(_field(m, "JP700_ma")) / _array(_field(m, "JP700_mj"))) * 100), "-b")
    ax.plot(
        x,
        _series(
            x,
            (
                1
                - _array(_field(m, "JP700_sa"))
                / (
                    _array(_field(m, "JP700_sjj")) * (_array(_field(m, "which_JP700_ma")) == 1)
                    + _array(_field(m, "JP700_scj")) * (_array(_field(m, "which_JP700_ma")) == 2)
                )
            )
            * 100,
        ),
        "-r",
    )
    ax.set_xlim(np.min(x), np.max(x))
    ax.set_ylim(0, 100)
    ax.set_ylabel("Source downregulation (%)")
    _annotate(ax, xpos, ypos, "(l)")

    ax = figure1.add_subplot(5, 4, 13)
    ax.plot(x, _series(x, _array(_field(m, "An_ma")) * 1e6), "-b")
    ax.plot(x, _series(x, _array(_field(m, "An_sa")) * 1e6), "-r")
    ax.set_xlim(np.min(x), np.max(x))
    ax.set_ylim(0, 30)
    ax.set_ylabel("An (umol CO2 m-2 s-1)")
    _annotate(ax, xpos, ypos, "(m)")

    ax = figure1.add_subplot(5, 4, 14)
    ax.plot(x, _series(x, _array(_field(m, "C_m")) / _array(_field(m, "O_m")) * 1000), "-b")
    ax.plot(x, _series(x, _array(_field(m, "C_sa")) / _array(_field(m, "O_sa")) * 1000), "-r")
    ax.set_xlim(np.min(x), np.max(x))
    ax.set_ylim(0, 30)
    ax.set_ylabel("CO2:O2 (mbar bar-1)")
    _annotate(ax, xpos, ypos, "(n)")

    ax = figure1.add_subplot(5, 4, 15)
    ax.plot(x, _series(x, _array(_field(m, "JP700_ma")) / _array(_field(m, "JP700_mc")) * _array(_field(m, "RUB_m")) * 1e6), "-b")
    ax.plot(
        x,
        _series(
            x,
            _array(_field(m, "JP700_sa"))
            / (
                _array(_field(m, "JP700_sjc")) * (_array(_field(m, "which_JP700_ma")) == 1)
                + _array(_field(m, "JP700_scc")) * (_array(_field(m, "which_JP700_ma")) == 2)
            )
            * _array(_field(m, "RUB_s"))
            * 1e6,
        ),
        "-r",
    )
    ax.set_xlim(np.min(x), np.max(x))
    ax.set_ylim(0, 15)
    ax.set_ylabel("Rub active sites (umol m-2)")
    _annotate(ax, xpos, ypos, "(o)")

    ax = figure1.add_subplot(5, 4, 16)
    ax.plot(x, _series(x, (1 - _array(_field(m, "JP700_ma")) / _array(_field(m, "JP700_mc"))) * 100), "-b")
    ax.plot(
        x,
        _series(
            x,
            (
                1
                - _array(_field(m, "JP700_sa"))
                / (
                    _array(_field(m, "JP700_sjc")) * (_array(_field(m, "which_JP700_ma")) == 1)
                    + _array(_field(m, "JP700_scc")) * (_array(_field(m, "which_JP700_ma")) == 2)
                )
            )
            * 100,
        ),
        "-r",
    )
    ax.set_xlim(np.min(x), np.max(x))
    ax.set_ylim(0, 100)
    ax.set_ylabel("Sink downregulation (%)")
    _annotate(ax, xpos, ypos, "(p)")

    ax = figure1.add_subplot(5, 4, 17)
    ax.plot(x, _series(x, _array(_field(m, "which_JP700_ma"))), "-b")
    ax.plot(
        x,
        _series(
            x,
            _array(_field(m, "which_JP700_sj")) * (_array(_field(m, "which_JP700_ma")) == 1)
            + _array(_field(m, "which_JP700_sc")) * (_array(_field(m, "which_JP700_ma")) == 2),
        ),
        "-r",
    )
    ax.set_xlim(np.min(x), np.max(x))
    ax.set_ylim(0, 3)
    ax.set_xlabel(xlab)
    ax.set_ylabel("Limiting states")
    xlim_curr = ax.get_xlim()
    ax.text(
        xlim_curr[0] + (xlim_curr[1] - xlim_curr[0]) * 0.5,
        0.66,
        "1: Light-limited",
        horizontalalignment="center",
    )
    ax.text(
        xlim_curr[0] + (xlim_curr[1] - xlim_curr[0]) * 0.5,
        0.33,
        "2: Light-saturated",
        horizontalalignment="center",
    )
    _annotate(ax, xpos, ypos, "(q)")

    ax = figure1.add_subplot(5, 4, 18)
    ax.plot(x, _series(x, _array(_field(m, "Kn2_ma")) / 1e9), "-b")
    ax.plot(x, _series(x, _array(_field(m, "Kn2_sa")) / 1e9), "-r")
    ax.set_xlim(np.min(x), np.max(x))
    ax.set_ylim(0, 6)
    ax.set_xlabel(xlab)
    ax.set_ylabel("Kn2 (ns-1)")
    _annotate(ax, xpos, ypos, "(r)")

    ax = figure1.add_subplot(5, 4, 19)
    ax.plot(
        x,
        _series(
            x,
            (_array(_field(m, "JP700_ma")) - _array(_field(m, "JP680_ma")))
            / _array(_field(m, "JP700_ma"))
            * 100,
        ),
        "-b",
    )
    ax.plot(
        x,
        _series(
            x,
            (_array(_field(m, "JP700_sa")) - _array(_field(m, "JP680_sa")))
            / _array(_field(m, "JP700_sa"))
            * 100,
        ),
        "-r",
    )
    ax.set_xlim(np.min(x), np.max(x))
    ax.set_ylim(0, 100)
    ax.set_xlabel(xlab)
    ax.set_ylabel("CEF1 fraction (%)")
    _annotate(ax, xpos, ypos, "(s)")

    ax = figure1.add_subplot(5, 4, 20)
    if _strcmp(_field(v, "pathway_opt"), "C3"):
        ax.plot(x, np.zeros(np.size(x)), "-r")
    else:
        ax.plot(
            x,
            _series(
                x,
                _array(_field(m, "L_C_sa")) / (_array(_field(m, "Vp_ma")) + _array(_field(m, "Vg_ma"))) * 100,
            ),
            "-r",
        )
    ax.set_xlim(np.min(x), np.max(x))
    ax.set_ylim(0, 100)
    ax.set_xlabel(xlab)
    ax.set_ylabel("L/(Vp+Vg) (%)")
    _annotate(ax, xpos, ypos, "(t)")

    return SimpleNamespace(figure1=figure1)
