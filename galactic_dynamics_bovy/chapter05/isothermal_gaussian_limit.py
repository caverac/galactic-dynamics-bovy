"""Convergence of a finite isotropic system to a Maxwellian velocity distribution.

Problem 5.9 shows analytically that an infinite isotropic equilibrium with a
spatially constant velocity dispersion must have a Gaussian (Maxwellian)
velocity distribution. This module demonstrates the result numerically with
galpy.

A King model is a finite isotropic system with

    f(E) ~ exp(E / sigma^2) - 1,   for E > 0,   f = 0 otherwise,

where E is the relative (binding) energy. It is finite precisely because f
vanishes for unbound energies, which truncates the local velocity distribution
at the escape speed v_esc(r) = sqrt(2 Psi(r)). The dimensionless central depth
W0 = Psi(0) / sigma^2 sets how far that truncation sits beyond the thermal core:
at the center v_esc / sigma = sqrt(2 W0). As W0 grows the King model approaches
the (infinite) isothermal sphere and its central velocity distribution
approaches a Gaussian.

For each W0 we evaluate the central velocity distribution exactly from the King
DF (galpy's ``kingdf.fE``) and quantify its shape with the excess kurtosis of a
single velocity component, which is zero for a Gaussian. The kurtosis is
computed from the speed moments using the isotropy identity

    kappa = (9/5) <v^4> / <v^2>^2 - 3,

since for an isotropic distribution <v_x^2> = <v^2>/3 and <v_x^4> = <v^4>/5.
"""

from pathlib import Path

from galpy.df import kingdf
from galpy.potential import evaluatePotentials, KingPotential
from matplotlib import pyplot as plt
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.ticker import AutoMinorLocator
import numpy as np

from galactic_dynamics_bovy.utils.assets import register_asset, save_figure_if_changed

# King concentrations from a barely-bound system to a near-isothermal one.
# galpy's King solver becomes unstable beyond W0 ~ 30, so we stay well below it.
W0_VALUES: tuple[float, ...] = (1.0, 2.0, 3.0, 4.0, 6.0, 8.0, 10.0, 13.0, 16.0, 20.0)


def central_excess_kurtosis(W0: float, *, rt: float = 5.0, n_grid: int = 4000) -> float:
    """Exact excess kurtosis of a King model's central velocity distribution.

    The central speed distribution is ``p(v) ~ v^2 f(Phi(0) + v^2/2)`` for
    ``v`` up to the escape speed, evaluated from galpy's King DF. The kurtosis
    of a single velocity component follows from the isotropy identity
    ``kappa = (9/5) <v^4>/<v^2>^2 - 3``.

    Parameters
    ----------
    W0 : float
        Dimensionless central potential depth Psi(0) / sigma^2.
    rt : float, optional
        Tidal (truncation) radius in natural units. Default is 5.0.
    n_grid : int, optional
        Number of speed grid points for the quadrature. Default is 4000.

    Returns
    -------
    float
        Excess kurtosis of a central velocity component (zero for a Gaussian).
    """
    kdf = kingdf(W0=W0, M=1.0, rt=rt)
    sigma = float(kdf.sigma)
    phi0 = float(evaluatePotentials(KingPotential(W0=W0, M=1.0, rt=rt), 1e-7, 0.0))

    v_esc = sigma * np.sqrt(2.0 * W0)
    v = np.linspace(0.0, v_esc, n_grid)
    f_e = np.clip(np.asarray(kdf.fE(phi0 + 0.5 * v**2), dtype=np.float64), 0.0, None)

    norm = float(np.trapezoid(v**2 * f_e, v))
    v2 = float(np.trapezoid(v**4 * f_e, v)) / norm
    v4 = float(np.trapezoid(v**6 * f_e, v)) / norm
    return 1.8 * v4 / v2**2 - 3.0


def _setup_axes(ax: Axes) -> None:
    """Apply standard tick and grid formatting to an axes object."""
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.tick_params(which="minor", length=3, color="gray", direction="in")
    ax.tick_params(which="major", length=6, direction="in")
    ax.tick_params(top=True, right=True, which="both")
    ax.grid(True, alpha=0.3)


@register_asset("p05_09_king_maxwellian_convergence.png")
def plot_king_maxwellian_convergence(path: Path | None = None) -> None:
    """Plot the central excess kurtosis versus W0.

    As W0 grows the King model approaches the isothermal sphere and the central
    velocity distribution approaches a Gaussian, so the excess kurtosis rises
    smoothly toward zero.

    Parameters
    ----------
    path : Path or None
        If provided, save the figure to this path. Otherwise display it.
    """
    w0 = np.array(W0_VALUES)
    kappa = np.array([central_excess_kurtosis(value) for value in W0_VALUES])

    fig: Figure
    ax: Axes
    fig, ax = plt.subplots(figsize=(6.0, 4.2))
    fig.subplots_adjust(left=0.14, right=0.96, top=0.94, bottom=0.14)

    ax.axhline(0.0, linestyle="--", color="#737373", linewidth=0.9)
    ax.plot(w0, kappa, "o-", color="#AA3377", linewidth=1.2, markersize=5)
    ax.set_xlabel(r"$W_0 = \Psi(0)\,/\,\sigma^2$")
    ax.set_ylabel(r"central excess kurtosis $\kappa$")
    ax.text(
        0.96,
        0.12,
        "Gaussian (isothermal) limit",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=9,
        color="#737373",
    )
    _setup_axes(ax)

    if path:
        save_figure_if_changed(fig, path, dpi=150, bbox_inches="tight")
        plt.close(fig)
    else:
        plt.show()


if __name__ == "__main__":  # pragma: no cover
    plot_king_maxwellian_convergence()
