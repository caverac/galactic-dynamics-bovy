"""Cored (non-singular) isothermal sphere density profile.

Problem 5.12 reduces the structure of a cored isothermal sphere to the
dimensionless, parameter-free equation for ``f(x) = rho(r)/rho_0`` with
``x = r/r_0`` and ``r_0^2 = 9 sigma^2 / (4 pi G rho_0)``. Writing
``y = ln f`` it reads

    (1/x^2) d/dx ( x^2 dy/dx ) = -9 e^y,   y(0) = 0,  y'(0) = 0.

This module integrates that equation outward from the core and compares the
resulting profile with the singular isothermal sphere ``rho = sigma^2 /
(2 pi G r^2)``, which in these units is ``f_sing(x) = 2 / (9 x^2)``.

The core boundary conditions fix the small-x series
``y = -(3/2) x^2 + (27/40) x^4 + ...``, used to start the integration just off
the origin (where the ``2/x`` term is singular).
"""

from dataclasses import dataclass
from pathlib import Path

from matplotlib import pyplot as plt
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.ticker import AutoMinorLocator
import numpy as np
import numpy.typing as npt
from scipy.integrate import solve_ivp

from galactic_dynamics_bovy.utils.assets import register_asset, save_figure_if_changed


@dataclass
class IsothermalProfile:
    """Dimensionless cored isothermal sphere profile.

    Attributes
    ----------
    x : ndarray
        Scaled radius r / r_0.
    f : ndarray
        Scaled density rho / rho_0 = f(x).
    """

    x: npt.NDArray[np.float64]
    f: npt.NDArray[np.float64]


def _isothermal_rhs(x: float, state: npt.NDArray[np.float64]) -> list[float]:
    """Right-hand side of the isothermal sphere ODE as a first-order system.

    With ``state = [y, y']`` and ``y = ln(rho/rho_0)``,
    ``y'' = -9 e^y - (2/x) y'``.

    Parameters
    ----------
    x : float
        Scaled radius.
    state : ndarray
        Two-element array ``[y, y']``.

    Returns
    -------
    list of float
        Derivatives ``[y', y'']``.
    """
    y, y_prime = state
    return [y_prime, -9.0 * np.exp(y) - 2.0 * y_prime / x]


def cored_isothermal_profile(
    *,
    x_min: float = 1e-2,
    x_max: float = 1e2,
    n_points: int = 400,
    x0: float = 1e-3,
) -> IsothermalProfile:
    """Integrate the cored isothermal sphere ODE to get f(x) = rho/rho_0.

    The integration starts at ``x0`` using the core series expansion
    ``y = -(3/2) x^2 + (27/40) x^4`` and runs out to ``x_max``, sampling on a
    logarithmic grid.

    Parameters
    ----------
    x_min : float, optional
        Smallest sampled scaled radius (>= x0). Default is 1e-2.
    x_max : float, optional
        Largest sampled scaled radius. Default is 100.
    n_points : int, optional
        Number of logarithmically spaced sample points. Default is 400.
    x0 : float, optional
        Radius at which the series initial condition is applied. Default 1e-3.

    Returns
    -------
    IsothermalProfile
        The grid ``x`` and the density profile ``f``.
    """
    y0 = -1.5 * x0**2 + (27.0 / 40.0) * x0**4
    y_prime0 = -3.0 * x0 + (27.0 / 10.0) * x0**3

    x = np.logspace(np.log10(x_min), np.log10(x_max), n_points)
    sol = solve_ivp(
        _isothermal_rhs,
        (x0, x_max),
        np.array([y0, y_prime0]),
        t_eval=x,
        rtol=1e-10,
        atol=1e-12,
    )
    return IsothermalProfile(x=x, f=np.exp(sol.y[0]))


def singular_isothermal_profile(x: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
    """Singular isothermal sphere density in scaled units: f_sing = 2 / (9 x^2).

    Parameters
    ----------
    x : ndarray
        Scaled radius r / r_0.

    Returns
    -------
    ndarray
        Scaled density rho / rho_0 for the singular isothermal sphere.
    """
    return 2.0 / (9.0 * x**2)


def _setup_axes(ax: Axes) -> None:
    """Apply standard tick and grid formatting to an axes object."""
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.tick_params(which="minor", length=3, color="gray", direction="in")
    ax.tick_params(which="major", length=6, direction="in")
    ax.tick_params(top=True, right=True, which="both")


@register_asset("p05_12_cored_isothermal_sphere.png")
def plot_cored_isothermal_sphere(path: Path | None = None) -> None:
    """Plot the cored isothermal sphere profile against the singular one.

    Parameters
    ----------
    path : Path or None
        If provided, save the figure to this path. Otherwise display it.
    """
    profile = cored_isothermal_profile()
    f_singular = singular_isothermal_profile(profile.x)

    fig: Figure
    ax: Axes
    fig, ax = plt.subplots(figsize=(6.0, 4.6))
    fig.subplots_adjust(left=0.14, right=0.96, top=0.95, bottom=0.13)

    ax.loglog(profile.x, profile.f, "-", color="#000000", linewidth=2, label="cored isothermal sphere")
    ax.loglog(profile.x, f_singular, "--", color="#444444", linewidth=2, label=r"singular: $2/(9x^2)$")
    ax.set_xlabel(r"$r\,/\,r_0$")
    ax.set_ylabel(r"$\rho\,/\,\rho_0$")
    ax.set_ylim(1e-5, 2.0)
    ax.axhline(1.0, color="#4A4A4A", linestyle=":", linewidth=1)
    ax.axvline(1.0, color="#4A4A4A", linestyle=":", linewidth=1)
    ax.legend(frameon=False, fontsize=9)
    _setup_axes(ax)

    if path:
        save_figure_if_changed(fig, path, dpi=150, bbox_inches="tight")
        plt.close(fig)
    else:
        plt.show()


if __name__ == "__main__":  # pragma: no cover
    plot_cored_isothermal_sphere()
