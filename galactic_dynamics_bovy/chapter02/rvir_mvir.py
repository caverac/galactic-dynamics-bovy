"""Virial radius and mass for NFW halos as a function of overdensity.

This module computes how the virial radius and virial mass of an NFW halo
vary with the overdensity threshold delta_v, which is defined relative to the
critical density rho_crit.

For an NFW profile with central density rho0 and scale radius a, the mean
density within radius r is:

    rho(<r) = 3 rho0 f(c) / c^3

where c = r/a is the concentration and f(c) = ln(1+c) - c/(1+c) is the
NFW mass function. The overdensity is then delta = rho(<r) / rho_crit.
"""

from pathlib import Path

from matplotlib import pyplot as plt
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.ticker import AutoMinorLocator
import numpy as np
import numpy.typing as npt

from galactic_dynamics_bovy.utils.assets import register_asset, save_figure_if_changed
from galactic_dynamics_bovy.utils.units import GalacticUnits


def _nfw_f(c: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
    """Compute the NFW mass function f(c) = ln(1+c) - c/(1+c).

    Parameters
    ----------
    c : ndarray
        Concentration parameter (r/a).

    Returns
    -------
    ndarray
        Values of f(c).
    """
    return np.log(1 + c) - c / (1 + c)


def get_rvir_mvir_delta(
    rho0: float,
    a: float,
    rho_crit: float,
    /,
    *,
    rvir_range: tuple[float, float] = (80, 250),
    n_points: int = 500,
) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64], npt.NDArray[np.float64]]:
    """Compute virial quantities as a function of overdensity for an NFW halo.

    Given NFW parameters (rho0, a), sweeps over a range of virial radii and
    computes the corresponding overdensity delta and enclosed mass M_vir at each
    radius.

    The overdensity is computed as:
        delta = 3 (rho0/rho_crit) f(r_vir/a) / (r_vir/a)^3

    And the virial mass is:
        M_vir = (4*pi/3) r_vir^3 rho_crit delta

    Parameters
    ----------
    rho0 : float
        NFW central density parameter rho0 (same units as rho_crit).
    a : float
        NFW scale radius (length units, e.g., kpc).
    rho_crit : float
        Critical density of the universe (same units as rho0).
    rvir_range : tuple[float, float]
        Range of virial radii to sweep (same units as a).
    n_points : int
        Number of points to evaluate.

    Returns
    -------
    rvir : ndarray
        Virial radius values (same units as a).
    mvir : ndarray
        Virial mass at each radius (mass units from rho0 a^3).
    delta : ndarray
        Overdensity delta = rho(<r_vir) / rho_crit at each radius.
    """
    rvir = np.linspace(*rvir_range, n_points)
    c = rvir / a

    delta = 3 * (rho0 / rho_crit) * _nfw_f(c) / c**3
    mvir = 4 * np.pi * rvir**3 * rho_crit * delta / 3

    return rvir, mvir, delta


@register_asset("p02_02_nfw_rvir_mvir_delta.png")
def plot_rvir_mvir_delta(path: Path | None = None) -> None:
    """Plot virial mass and radius vs overdensity for an NFW halo.

    Creates a two-panel plot showing how M_vir and r_vir depend on the
    overdensity threshold delta_v for a fixed NFW density profile.

    Parameters
    ----------
    path : Path or None
        If provided, save the figure to this path. Otherwise display it.
    """
    units = GalacticUnits(h=1.0)

    # NFW parameters in GalacticUnits: (kpc, 10^10 Msun, Gyr)
    rho0 = 0.00035  # 10^10 Msun / kpc^3
    a = 16.0  # kpc

    rvir, mvir, delta = get_rvir_mvir_delta(rho0, a, units.rho_crit)

    fig: Figure
    axs: list[Axes]
    fig, axs = plt.subplots(2, 1, figsize=(6.5, 7), sharex=True)

    fig.subplots_adjust(
        left=0.1,
        right=0.92,
        bottom=0.1,
        top=0.92,
        wspace=0.05,
        hspace=0.05,
    )

    for ax in axs:
        ax.xaxis.set_minor_locator(AutoMinorLocator())
        ax.yaxis.set_minor_locator(AutoMinorLocator())
        ax.tick_params(which="minor", length=3, color="gray", direction="in")
        ax.tick_params(which="major", length=6, direction="in")
        ax.tick_params(top=True, right=True, which="both")

    # plot mvir vs delta
    ax = axs[0]
    ax.plot(delta, mvir, "-", color="#000000", linewidth=2)
    ax.set_ylabel(r"$M_{\mathrm{vir}}$ ($10^{10}M_\odot$)")
    ax.set_ylim(mvir.min() * 0.95, mvir.max() * 1.05)
    ax.grid(True, alpha=0.3)

    # plot rvir vs delta
    ax = axs[1]
    ax.plot(delta, rvir, "-", color="#000000", linewidth=2)
    ax.set_xlabel(r"$\Delta_v$")
    ax.set_ylabel(r"$r_{\mathrm{vir}}$ (kpc)")
    ax.set_ylim(rvir.min() * 0.95, rvir.max() * 1.05)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()

    if path:
        save_figure_if_changed(fig, path, dpi=150, bbox_inches="tight")
        plt.close(fig)
    else:
        plt.show()
