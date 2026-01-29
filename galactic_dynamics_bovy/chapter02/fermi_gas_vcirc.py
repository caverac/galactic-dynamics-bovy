"""Circular velocity for a self-gravitating degenerate Fermi gas.

This module computes the circular velocity curve for a self-gravitating
degenerate Fermi gas (polytrope with n=1), as derived in Problem 2.13.

The density profile is:
    ρ(r) = ρ_c sinc(kr)

where k = π/R and R is the boundary radius where ρ(R) = 0.

The circular velocity is:
    v_c² = (4πGρ_c/k²) * (sinc(kr) - cos(kr))  for r ≤ R

The total mass is:
    M = 4ρ_c R³/π

So given M and R, we can derive:
    ρ_c = πM/(4R³)
    k = π/R
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


def get_fermi_gas_params(
    total_mass: float,
    boundary_radius: float,
) -> tuple[float, float]:
    """Convert total mass and boundary radius to Fermi gas parameters.

    Parameters
    ----------
    total_mass : float
        Total mass M in units of 10^10 Msun.
    boundary_radius : float
        Boundary radius R in kpc.

    Returns
    -------
    rho_c : float
        Central density in 10^10 Msun / kpc³.
    k : float
        Wave number in 1/kpc.
    """
    k = np.pi / boundary_radius
    rho_c = np.pi * total_mass / (4 * boundary_radius**3)
    return rho_c, k


def vcirc_fermi_gas(
    r: npt.NDArray[np.float64],
    rho_c: float,
    k: float,
    G: float,
) -> npt.NDArray[np.float64]:
    """Compute circular velocity for a Fermi gas profile.

    Parameters
    ----------
    r : ndarray
        Radius in kpc.
    rho_c : float
        Central density in 10^10 Msun / kpc³.
    k : float
        Wave number in 1/kpc.
    G : float
        Gravitational constant in kpc (km/s)² / (10^10 Msun).

    Returns
    -------
    ndarray
        Circular velocity in km/s.
    """
    kr = k * r
    # Handle r=0 case: sinc(0) - cos(0) = 1 - 1 = 0, so v_c(0) = 0
    with np.errstate(divide="ignore", invalid="ignore"):
        sinc_kr = np.sinc(kr / np.pi)  # np.sinc(x) = sin(πx)/(πx)
        term = sinc_kr - np.cos(kr)

    # v_c² = (4πGρ_c/k²) * (sinc(kr) - cos(kr))
    vc_squared = 4 * np.pi * G * rho_c / k**2 * term
    vc_squared = np.maximum(vc_squared, 0)  # Ensure non-negative

    return np.sqrt(vc_squared)


@register_asset("fermi_gas_vcirc.png")
def plot_fermi_gas_vcirc(*, path: Path | None = None) -> None:
    """Plot circular velocity for a Fermi gas with M=1e11 Msun, R=15 kpc.

    Parameters
    ----------
    path : Path or None
        If provided, save the figure to this path. Otherwise display it.
    """
    units = GalacticUnits()

    # Parameters: M = 1e11 Msun = 10 (in 10^10 Msun), R = 15 kpc
    total_mass = 10.0  # 10^10 Msun
    boundary_radius = 15.0  # kpc

    rho_c, k = get_fermi_gas_params(total_mass, boundary_radius)

    # Radius grid from 0 to R
    r = np.linspace(0, boundary_radius, 500)
    vc = vcirc_fermi_gas(r, rho_c, k, units.G_kms)

    fig: Figure
    ax: Axes
    fig, ax = plt.subplots(figsize=(6.5, 4.5))

    ax.plot(r, vc, "-", color="#000000", linewidth=2)

    ax.set_xlabel(r"$r$ (kpc)")
    ax.set_ylabel(r"$v_c$ (km/s)")
    ax.set_xlim(0, boundary_radius)
    ax.set_ylim(0, None)

    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.tick_params(which="minor", length=3, color="gray", direction="in")
    ax.tick_params(which="major", length=6, direction="in")
    ax.tick_params(top=True, right=True, which="both")

    ax.grid(True, alpha=0.3)

    plt.tight_layout()

    if path:
        save_figure_if_changed(fig, path, dpi=150, bbox_inches="tight")
        plt.close(fig)
    else:
        plt.show()
