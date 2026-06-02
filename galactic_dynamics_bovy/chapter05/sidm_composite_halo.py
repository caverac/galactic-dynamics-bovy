"""Composite SIDM halo: cored isothermal core matched to an outer NFW profile.

Problem 5.12 (part e). A simple model of a strongly self-interacting dark-matter
halo replaces the cuspy NFW centre with the cored isothermal sphere of parts
a-b inside the once-scattered radius r_1 (beta = 0 case, Problem 5.12d), keeping
NFW outside. The cored profile

    rho(r) = rho_0 f(r/r_0),    M(<r) = 4 pi rho_0 r_0^3 mu(r/r_0),

with f the universal cored solution and mu(X) = int_0^X f(x) x^2 dx, has two free
parameters rho_0 and r_0. We fix them by matching both the density and the
enclosed mass to the NFW halo at r_1. Writing X_1 = r_1/r_0 this reduces to a
single root equation,

    mu(X_1) / (X_1^3 f(X_1)) = M_NFW(<r_1) / (4 pi rho_NFW(r_1) r_1^3),

after which r_0 = r_1/X_1 and rho_0 = rho_NFW(r_1)/f(X_1).
"""

from dataclasses import dataclass
from pathlib import Path

from matplotlib import pyplot as plt
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.ticker import AutoMinorLocator
import numpy as np
import numpy.typing as npt
from scipy.integrate import cumulative_trapezoid
from scipy.optimize import brentq

from galactic_dynamics_bovy.chapter05.cored_isothermal_sphere import cored_isothermal_profile
from galactic_dynamics_bovy.chapter05.nfw_jeans_dispersion import NFWHalo
from galactic_dynamics_bovy.chapter05.sidm_scatter_radius import scatter_radius
from galactic_dynamics_bovy.utils.assets import register_asset, save_figure_if_changed


@dataclass
class CoredMatch:
    """Cored-isothermal parameters matched to an NFW halo at r_1.

    Attributes
    ----------
    rho_0 : float
        Central density of the cored isothermal sphere in Msun/kpc^3.
    r_0 : float
        Core scale radius in kpc.
    r_1 : float
        Matching (once-scattered) radius in kpc.
    x_match : float
        Dimensionless matching radius X_1 = r_1 / r_0.
    """

    rho_0: float
    r_0: float
    r_1: float
    x_match: float


def _universal_cored_solution(
    n_points: int = 4000,
) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64], npt.NDArray[np.float64]]:
    """Return (x, f, mu) for the universal cored isothermal sphere.

    ``f`` is the dimensionless density and ``mu(x) = int_0^x f x'^2 dx'`` the
    dimensionless enclosed mass.

    Parameters
    ----------
    n_points : int, optional
        Number of radius samples. Default is 4000.

    Returns
    -------
    tuple of ndarray
        The grids ``x``, ``f``, and ``mu``.
    """
    profile = cored_isothermal_profile(x_min=1e-3, x_max=1e3, n_points=n_points)
    integral = cumulative_trapezoid(profile.f * profile.x**2, profile.x, initial=0)
    mu = np.asarray(integral + profile.x[0] ** 3 / 3.0, dtype=np.float64)
    return profile.x, profile.f, mu


def match_cored_isothermal_to_nfw(halo: NFWHalo, r_1: float) -> CoredMatch:
    """Match a cored isothermal sphere to an NFW halo in density and mass at r_1.

    Parameters
    ----------
    halo : NFWHalo
        The outer NFW halo.
    r_1 : float
        Matching radius in kpc.

    Returns
    -------
    CoredMatch
        The fitted central density, core radius, and matching radius.
    """
    x, f, mu = _universal_cored_solution()
    log_x = np.log(x)

    def f_at(scaled_radius: float) -> float:
        return float(np.exp(np.interp(np.log(scaled_radius), log_x, np.log(f))))

    def mu_at(scaled_radius: float) -> float:
        return float(np.exp(np.interp(np.log(scaled_radius), log_x, np.log(mu))))

    rho_1 = halo.density(r_1)
    mass_1 = halo.enclosed_mass(r_1)
    target = mass_1 / (4.0 * np.pi * rho_1 * r_1**3)

    def excess(scaled_radius: float) -> float:
        return mu_at(scaled_radius) / (scaled_radius**3 * f_at(scaled_radius)) - target

    x_match = float(brentq(excess, 1e-2, 10.0))
    r_0 = r_1 / x_match
    rho_0 = rho_1 / f_at(x_match)
    return CoredMatch(rho_0=rho_0, r_0=r_0, r_1=r_1, x_match=x_match)


def composite_density(halo: NFWHalo, match: CoredMatch, radii: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
    """Composite SIDM density: cored isothermal inside r_1, NFW outside.

    Parameters
    ----------
    halo : NFWHalo
        The outer NFW halo.
    match : CoredMatch
        The matched cored-isothermal parameters.
    radii : ndarray
        Radii in kpc.

    Returns
    -------
    ndarray
        Density in Msun/kpc^3.
    """
    x, f, _ = _universal_cored_solution()
    log_x = np.log(x)
    cored = match.rho_0 * np.exp(np.interp(np.log(radii / match.r_0), log_x, np.log(f)))
    nfw = np.array([halo.density(float(r)) for r in radii])
    return np.where(radii <= match.r_1, cored, nfw)


def _setup_axes(ax: Axes) -> None:
    """Apply standard tick formatting to an axes object."""
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.tick_params(which="minor", length=3, color="gray", direction="in")
    ax.tick_params(which="major", length=6, direction="in")
    ax.tick_params(top=True, right=True, which="both")


@register_asset("p05_12_sidm_composite_halo.png")
def plot_sidm_composite_halo(path: Path | None = None) -> None:
    """Plot the composite SIDM density profile against the pure NFW profile.

    Parameters
    ----------
    path : Path or None
        If provided, save the figure to this path. Otherwise display it.
    """
    halo = NFWHalo.from_virial_mass()
    r_1 = scatter_radius(halo, 0.0)
    match = match_cored_isothermal_to_nfw(halo, r_1)

    radii = np.logspace(0.0, np.log10(300.0), 300)
    composite = composite_density(halo, match, radii)
    nfw = np.array([halo.density(float(r)) for r in radii])

    fig: Figure
    ax: Axes
    fig, ax = plt.subplots(figsize=(6.0, 4.6))
    fig.subplots_adjust(left=0.14, right=0.96, top=0.95, bottom=0.13)

    ax.loglog(radii, nfw, "--", color="#888888", linewidth=1.4, label="NFW (cuspy)")
    ax.loglog(radii, composite, "-", color="#000000", linewidth=2, label="SIDM (cored + NFW)")
    ax.axvline(r_1, color="#737373", linestyle=":", linewidth=1.0)
    ax.text(r_1, ax.get_ylim()[1], r"$r_1$", color="#737373", ha="center", va="bottom")
    ax.set_xlabel(r"$r\; (\mathrm{kpc})$")
    ax.set_ylabel(r"$\rho\; (M_\odot\,\mathrm{kpc^{-3}})$")
    ax.set_xlim(1.0, 300.0)
    ax.legend(frameon=False)
    _setup_axes(ax)

    if path:
        save_figure_if_changed(fig, path, dpi=150, bbox_inches="tight")
        plt.close(fig)
    else:
        plt.show()


if __name__ == "__main__":  # pragma: no cover
    plot_sidm_composite_halo()
