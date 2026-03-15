from types import SimpleNamespace

import matplotlib.pyplot as plt
import numpy as np


def _field(value, name):
    if isinstance(value, dict):
        return value[name]
    return getattr(value, name)


def _strcmp(pathway_option, target):
    if isinstance(pathway_option, (list, tuple)) and len(pathway_option) == 1:
        pathway_option = pathway_option[0]
    return pathway_option == target


def _array(value):
    return np.asarray(value, dtype=float)


def _series(x, value):
    value = np.asarray(value, dtype=float)
    if value.ndim == 0 or value.size == 1:
        return np.full(np.asarray(x).shape, float(value))
    return value


def _safe_divide(numerator, denominator, fill_value=0.0):
    numerator = np.asarray(numerator, dtype=float)
    denominator = np.asarray(denominator, dtype=float)
    numerator, denominator = np.broadcast_arrays(numerator, denominator)
    result = np.full(numerator.shape, fill_value, dtype=float)
    np.divide(numerator, denominator, out=result, where=denominator != 0)
    return result


def _annotate(ax, xpos, ypos, label):
    xlim_curr = ax.get_xlim()
    ylim_curr = ax.get_ylim()
    ax.text(
        xlim_curr[0] + (xlim_curr[1] - xlim_curr[0]) * xpos,
        ylim_curr[0] + (ylim_curr[1] - ylim_curr[0]) * ypos,
        label,
    )


def plotter_forward_fun_c3c4(outputname, v, m):
    if np.ptp(_array(_field(m, "PPFD"))) > 0:
        x = _array(_field(m, "PPFD")) * 1e6
        xlab = "PPFD (umol PAR m-2 s-1)"

    if np.ptp(_array(_field(m, "CO2_m"))) > 0:
        x = _array(_field(m, "CO2_m")) * 1e6
        xlab = "CO2_m (ubar CO2)"

    if np.ptp(_array(_field(m, "Temp"))) > 0:
        x = _array(_field(m, "Temp"))
        xlab = "Temp (degrees C)"

    if np.ptp(_array(_field(m, "PPFD"))) > 0 and np.ptp(_array(_field(m, "Temp"))) > 0:
        x = np.linspace(1, 1440, 1440) / 60
        xlab = "Time (hours)"

    figure1 = plt.figure()
    xpos = 0.85
    ypos = 0.9

    ax = figure1.add_subplot(5, 4, 1)
    ax.plot(x, _series(x, _array(_field(m, "J_PSII_m_actual")) * 1e6), "-b")
    ax.plot(x, _series(x, _array(_field(m, "J_PSII_s_actual")) * 1e6), "-r")
    ax.set_xlim(np.min(x), np.max(x))
    ax.set_ylim(0, 250)
    ax.set_ylabel("J_PSII_*_actual (umol e- m-2 s-1)")
    ax.legend(["Mesophyll", "Bundle sheath"])
    _annotate(ax, xpos, ypos, "(a)")

    ax = figure1.add_subplot(5, 4, 2)
    ax.plot(x, _series(x, _array(_field(m, "phi_P2_m_actual")) * 100), "-b")
    ax.plot(x, _series(x, _array(_field(m, "phi_P2_s_actual")) * 100), "-r")
    ax.set_xlim(np.min(x), np.max(x))
    ax.set_ylim(0, 100)
    ax.set_ylabel("phi_P2_*_actual (%)")
    _annotate(ax, xpos, ypos, "(b)")

    ax = figure1.add_subplot(5, 4, 3)
    ax.plot(x, _series(x, _array(_field(m, "phi_N2_m_actual")) * 100), "-b")
    ax.plot(x, _series(x, _array(_field(m, "phi_N2_s_actual")) * 100), "-r")
    ax.set_xlim(np.min(x), np.max(x))
    ax.set_ylim(0, 100)
    ax.set_ylabel("phi_N2_*_actual (%)")
    _annotate(ax, xpos, ypos, "(c)")

    ax = figure1.add_subplot(5, 4, 4)
    ax.plot(x, _series(x, (_array(_field(m, "phi_D2_m_actual")) + _array(_field(m, "phi_F2_m_actual"))) * 100), "-b")
    ax.plot(x, _series(x, (_array(_field(m, "phi_D2_s_actual")) + _array(_field(m, "phi_F2_s_actual"))) * 100), "-r")
    ax.set_xlim(np.min(x), np.max(x))
    ax.set_ylim(0, 100)
    ax.set_ylabel("phi_D2_*_actual + phi_F2_*_actual (%)")
    _annotate(ax, xpos, ypos, "(d)")

    ax = figure1.add_subplot(5, 4, 5)
    ax.plot(x, _series(x, _array(_field(m, "J_PSI_m_actual")) * 1e6), "-b")
    ax.plot(x, _series(x, _array(_field(m, "J_PSI_s_actual")) * 1e6), "-r")
    ax.set_xlim(np.min(x), np.max(x))
    ax.set_ylim(0, 250)
    ax.set_ylabel("J_PSI_*_actual (umol e- m-2 s-1)")
    _annotate(ax, xpos, ypos, "(e)")

    ax = figure1.add_subplot(5, 4, 6)
    ax.plot(x, _series(x, _array(_field(m, "phi_P1_m_actual")) * 100), "-b")
    ax.plot(x, _series(x, _array(_field(m, "phi_P1_s_actual")) * 100), "-r")
    ax.set_xlim(np.min(x), np.max(x))
    ax.set_ylim(0, 100)
    ax.set_ylabel("phi_P1_*_actual (%)")
    _annotate(ax, xpos, ypos, "(f)")

    ax = figure1.add_subplot(5, 4, 7)
    ax.plot(x, _series(x, _array(_field(m, "phi_N1_m_actual")) * 100), "-b")
    ax.plot(x, _series(x, _array(_field(m, "phi_N1_s_actual")) * 100), "-r")
    ax.set_xlim(np.min(x), np.max(x))
    ax.set_ylim(0, 100)
    ax.set_ylabel("phi_N1_*_actual (%)")
    _annotate(ax, xpos, ypos, "(g)")

    ax = figure1.add_subplot(5, 4, 8)
    ax.plot(x, _series(x, (_array(_field(m, "phi_D1_m_actual")) + _array(_field(m, "phi_F1_m_actual"))) * 100), "-b")
    ax.plot(x, _series(x, (_array(_field(m, "phi_D1_s_actual")) + _array(_field(m, "phi_F1_s_actual"))) * 100), "-r")
    ax.set_xlim(np.min(x), np.max(x))
    ax.set_ylim(0, 100)
    ax.set_ylabel("phi_D1_*_actual + phi_F1_*_actual (%)")
    _annotate(ax, xpos, ypos, "(h)")

    ax = figure1.add_subplot(5, 4, 9)
    ax.plot(x, _series(x, _array(_field(m, "J_PSI_m_actual")) * 1e6), "-b")
    ax.plot(x, _series(x, _array(_field(m, "J_PSI_s_actual")) * 1e6), "-r")
    ax.set_xlim(np.min(x), np.max(x))
    ax.set_ylim(0, 250)
    ax.set_ylabel("Cyt b6f ETR (umol e- m-2 s-1)")
    _annotate(ax, xpos, ypos, "(i)")

    ax = figure1.add_subplot(5, 4, 10)
    ax.plot(x, _series(x, (1 - _array(_field(m, "q_P2_m_actual"))) * 100), "-b")
    ax.plot(x, _series(x, (1 - _array(_field(m, "q_P2_s_actual"))) * 100), "-r")
    ax.set_xlim(np.min(x), np.max(x))
    ax.set_ylim(0, 100)
    ax.set_ylabel("1 - q_P2_*_actual (%)")
    _annotate(ax, xpos, ypos, "(j)")

    ax = figure1.add_subplot(5, 4, 11)
    ax.plot(
        x,
        _series(
            x,
            _safe_divide(_array(_field(m, "J_PSI_m_actual")), _array(_field(m, "J_PSI_m_j")))
            * _safe_divide(_array(_field(m, "V_q_max_m")), _array(_field(m, "Cytbf_density_m"))),
        ),
        "-b",
    )
    ax.plot(
        x,
        _series(
            x,
            _safe_divide(
                _array(_field(m, "J_PSI_s_actual")),
                _array(_field(m, "J_PSI_s_jj")) * (_array(_field(m, "which_J_PSI_m")) == 1)
                + _array(_field(m, "J_PSI_s_cj")) * (_array(_field(m, "which_J_PSI_m")) == 2),
            )
            * _safe_divide(_array(_field(m, "V_q_max_s")), _array(_field(m, "Cytbf_density_s"))),
        ),
        "-r",
    )
    ax.set_xlim(np.min(x), np.max(x))
    ax.set_ylim(0, 600)
    ax.set_ylabel("Cyt b6f rate constant (s-1)")
    _annotate(ax, xpos, ypos, "(k)")

    ax = figure1.add_subplot(5, 4, 12)
    ax.plot(x, _series(x, (1 - _array(_field(m, "J_PSI_m_actual")) / _array(_field(m, "J_PSI_m_j"))) * 100), "-b")
    ax.plot(
        x,
        _series(
            x,
            (
                1
                - _safe_divide(
                    _array(_field(m, "J_PSI_s_actual")),
                    _array(_field(m, "J_PSI_s_jj")) * (_array(_field(m, "which_J_PSI_m")) == 1)
                    + _array(_field(m, "J_PSI_s_cj")) * (_array(_field(m, "which_J_PSI_m")) == 2),
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
    ax.plot(x, _series(x, _array(_field(m, "A_net_m_actual")) * 1e6), "-b")
    ax.plot(x, _series(x, _array(_field(m, "A_net_s_actual")) * 1e6), "-r")
    ax.set_xlim(np.min(x), np.max(x))
    ax.set_ylim(0, 30)
    ax.set_ylabel("A_net_*_actual (umol CO2 m-2 s-1)")
    _annotate(ax, xpos, ypos, "(m)")

    ax = figure1.add_subplot(5, 4, 14)
    ax.plot(x, _series(x, _safe_divide(_array(_field(m, "CO2_m")), _array(_field(m, "O2_m"))) * 1000), "-b")
    ax.plot(x, _series(x, _safe_divide(_array(_field(m, "CO2_s_actual")), _array(_field(m, "O2_s_actual"))) * 1000), "-r")
    ax.set_xlim(np.min(x), np.max(x))
    ax.set_ylim(0, 30)
    ax.set_ylabel("CO2_* / O2_* (mbar bar-1)")
    _annotate(ax, xpos, ypos, "(n)")

    ax = figure1.add_subplot(5, 4, 15)
    ax.plot(
        x,
        _series(
            x,
            _safe_divide(_array(_field(m, "J_PSI_m_actual")), _array(_field(m, "J_PSI_m_c")))
            * _array(_field(m, "Rubisco_density_m"))
            * 1e6,
        ),
        "-b",
    )
    ax.plot(
        x,
        _series(
            x,
            _safe_divide(
                _array(_field(m, "J_PSI_s_actual")),
                _array(_field(m, "J_PSI_s_jc")) * (_array(_field(m, "which_J_PSI_m")) == 1)
                + _array(_field(m, "J_PSI_s_cc")) * (_array(_field(m, "which_J_PSI_m")) == 2),
            )
            * _array(_field(m, "Rubisco_density_s"))
            * 1e6,
        ),
        "-r",
    )
    ax.set_xlim(np.min(x), np.max(x))
    ax.set_ylim(0, 15)
    ax.set_ylabel("Rub active sites (umol m-2)")
    _annotate(ax, xpos, ypos, "(o)")

    ax = figure1.add_subplot(5, 4, 16)
    ax.plot(x, _series(x, (1 - _array(_field(m, "J_PSI_m_actual")) / _array(_field(m, "J_PSI_m_c"))) * 100), "-b")
    ax.plot(
        x,
        _series(
            x,
            (
                1
                - _safe_divide(
                    _array(_field(m, "J_PSI_s_actual")),
                    _array(_field(m, "J_PSI_s_jc")) * (_array(_field(m, "which_J_PSI_m")) == 1)
                    + _array(_field(m, "J_PSI_s_cc")) * (_array(_field(m, "which_J_PSI_m")) == 2),
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
    ax.plot(x, _series(x, _array(_field(m, "which_J_PSI_m"))), "-b")
    ax.plot(
        x,
        _series(
            x,
            _array(_field(m, "which_J_PSI_s_j")) * (_array(_field(m, "which_J_PSI_m")) == 1)
            + _array(_field(m, "which_J_PSI_s_c")) * (_array(_field(m, "which_J_PSI_m")) == 2),
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
    ax.plot(x, _series(x, _array(_field(m, "k_N2_m_actual")) / 1e9), "-b")
    ax.plot(x, _series(x, _array(_field(m, "k_N2_s_actual")) / 1e9), "-r")
    ax.set_xlim(np.min(x), np.max(x))
    ax.set_ylim(0, 6)
    ax.set_xlabel(xlab)
    ax.set_ylabel("k_N2_*_actual (ns-1)")
    _annotate(ax, xpos, ypos, "(r)")

    ax = figure1.add_subplot(5, 4, 19)
    ax.plot(
        x,
        _series(
            x,
            _safe_divide(
                _array(_field(m, "J_PSI_m_actual")) - _array(_field(m, "J_PSII_m_actual")),
                _array(_field(m, "J_PSI_m_actual")),
            )
            * 100,
        ),
        "-b",
    )
    ax.plot(
        x,
        _series(
            x,
            _safe_divide(
                _array(_field(m, "J_PSI_s_actual")) - _array(_field(m, "J_PSII_s_actual")),
                _array(_field(m, "J_PSI_s_actual")),
            )
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
    if _strcmp(_field(v, "pathway_option"), "C3"):
        ax.plot(x, np.zeros(np.size(x)), "-r")
    else:
        ax.plot(
            x,
            _series(
                x,
                _array(_field(m, "Leak_CO2_s_actual")) / (_array(_field(m, "V_p_m_actual")) + _array(_field(m, "V_g_m_actual"))) * 100,
            ),
            "-r",
        )
    ax.set_xlim(np.min(x), np.max(x))
    ax.set_ylim(0, 100)
    ax.set_xlabel(xlab)
    ax.set_ylabel("L/(Vp+Vg) (%)")
    _annotate(ax, xpos, ypos, "(t)")

    return SimpleNamespace(figure1=figure1)
