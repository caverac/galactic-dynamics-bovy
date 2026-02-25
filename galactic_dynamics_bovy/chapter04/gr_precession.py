"""GR apsidal precession for Keplerian orbits with the Schwarzschild correction.

Computes the apsidal precession per orbit arising from the first-order
GR correction to the Newtonian effective potential:

    Phi_eff(r) = -GM/r + L^2/(2r^2) - GM L^2/(c^2 r^3)

The precession is evaluated by numerical quadrature of the azimuthal-advance
integral with the regularising change of variables from Problem 4.9:

    2r = (r_p + r_a) + (r_a - r_p) sin(t),    t in [-pi/2, pi/2].
"""

from pathlib import Path

from matplotlib import pyplot as plt
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.ticker import AutoMinorLocator
import numpy as np
import numpy.typing as npt
from rich.console import Console
from scipy import integrate

from galactic_dynamics_bovy.utils.assets import register_asset, save_figure_if_changed
from galactic_dynamics_bovy.utils.units import SolarUnits


def _solve_EL(
    r_p: float,
    r_a: float,
    GM: float,
    c: float,
) -> tuple[float, float]:
    """Solve for orbital energy and squared angular momentum from turning points.

    Given pericenter *r_p* and apocenter *r_a*, finds *E* and *L^2* such that
    ``E = Phi_eff(r_p) = Phi_eff(r_a)``.

    Parameters
    ----------
    r_p : float
        Pericenter distance (m).
    r_a : float
        Apocenter distance (m).
    GM : float
        Gravitational parameter (m^3/s^2).
    c : float
        Speed of light (m/s).

    Returns
    -------
    E : float
        Orbital energy per unit mass (m^2/s^2).
    L2 : float
        Squared angular momentum per unit mass (m^4/s^2).
    """
    u_p = 1.0 / r_p
    u_a = 1.0 / r_a
    c2 = c * c

    # From E = Phi_eff(r_p) = Phi_eff(r_a), factor out (u_p - u_a):
    #   GM = L^2 [ (u_p + u_a)/2 - (GM/c^2)(u_p^2 + u_p u_a + u_a^2) ]
    denom = (u_p + u_a) / 2.0 - GM / c2 * (u_p**2 + u_p * u_a + u_a**2)
    L2 = GM / denom

    E = -GM * u_p + L2 * u_p**2 / 2.0 - GM * L2 * u_p**3 / c2
    return E, L2


def _integrand(
    t: float,
    S: float,
    D: float,
    LD: float,
    kinetic_coeffs: tuple[float, float, float, float],
) -> float:
    """Regularised integrand for the azimuthal-advance quadrature.

    The radial kinetic energy is expressed as a polynomial in u = 1/r:

        2(E - Phi_eff) = k0 + k1*u + k2*u^2 + k3*u^3

    where the coefficients are precomputed by the caller.
    """
    k0, k1, k2, k3 = kinetic_coeffs
    r = 0.5 * (S + D * np.sin(t))
    u = 1.0 / r
    u2 = u * u
    kinetic = k0 + u * (k1 + u * (k2 + k3 * u))
    if kinetic <= 0.0:
        return 0.0
    return float(LD * np.cos(t) * u2 / np.sqrt(kinetic))


def compute_delta_psi(
    r_p: float,
    e: float,
    /,
    *,
    GM: float = SolarUnits.GM,
    c: float = SolarUnits.c,
) -> float:
    """Compute the azimuthal advance per radial period.

    Uses the regularised substitution
    ``2r = (r_p + r_a) + (r_a - r_p) sin t``
    so that the integrand is smooth on ``[-pi/2, pi/2]``.

    Parameters
    ----------
    r_p : float
        Pericenter distance (m).
    e : float
        Orbital eccentricity (0 < e < 1).
    GM : float
        Gravitational parameter (m^3/s^2).
    c : float
        Speed of light (m/s).

    Returns
    -------
    float
        Azimuthal advance Delta-psi (rad).
    """
    r_a = r_p * (1.0 + e) / (1.0 - e)
    E, L2 = _solve_EL(r_p, r_a, GM, c)
    S = r_p + r_a
    D = r_a - r_p
    LD = np.sqrt(L2) * D
    c2 = c * c
    kinetic_coeffs = (2.0 * E, 2.0 * GM, -L2, 2.0 * GM * L2 / c2)

    result, _ = integrate.quad(_integrand, -np.pi / 2, np.pi / 2, args=(S, D, LD, kinetic_coeffs), limit=200)
    return result


def compute_precession_curve(
    a_array: npt.NDArray[np.float64],
    e: float,
    /,
    *,
    GM: float = SolarUnits.GM,
    c: float = SolarUnits.c,
) -> npt.NDArray[np.float64]:
    """Compute the precession (Delta-psi - 2pi) for an array of semi-major axes.

    Parameters
    ----------
    a_array : ndarray
        Semi-major axes (m).
    e : float
        Orbital eccentricity.
    GM : float
        Gravitational parameter (m^3/s^2).
    c : float
        Speed of light (m/s).

    Returns
    -------
    ndarray
        Precession Delta-psi - 2pi (rad) at each a.
    """
    return np.array([compute_delta_psi(a * (1.0 - e), e, GM=GM, c=c) - 2.0 * np.pi for a in a_array])


def _setup_axes(ax: Axes) -> None:
    """Apply standard formatting to an axes object."""
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.tick_params(which="minor", length=3, color="gray", direction="in")
    ax.tick_params(which="major", length=6, direction="in")
    ax.tick_params(top=True, right=True, which="both")
    ax.grid(True, alpha=0.3)


@register_asset("p04_09_gr_precession.png")
def plot_gr_precession(path: Path | None = None) -> None:
    """Plot GR apsidal precession vs pericenter distance for the solar system.

    Three eccentricities are shown; Mercury's (e = 0.205630) is highlighted.

    Parameters
    ----------
    path : Path or None
        If provided, save the figure to this path. Otherwise display it.
    """
    eccentricities: list[float] = [0.01, 0.205630, 0.50]
    labels: list[str] = [r"$e = 0.01$", r"$e = 0.2056$ (Mercury)", r"$e = 0.50$"]
    colors: list[str] = ["#4E4E4E", "#000000", "#434343"]
    styles: list[str] = ["--", "-", ":"]

    a_au = np.linspace(0.05, 1.0, 200)
    a_m = a_au * SolarUnits.AU

    fig: Figure
    ax: Axes
    fig, ax = plt.subplots(figsize=(6.0, 4.0))
    fig.subplots_adjust(hspace=0.08, wspace=0.08, left=0.10, right=0.98, top=0.98, bottom=0.1)

    for e, label, color, ls in zip(eccentricities, labels, colors, styles):
        # Numerical (thick)
        precession_rad = compute_precession_curve(a_m, e)
        precession_arcsec = np.degrees(precession_rad) * 3600.0
        ax.plot(a_au, precession_arcsec, ls, color=color, linewidth=2.5, label=label)

        # Analytical: 6*pi*GM / (c^2 * a * (1-e^2))
        analytical_rad = 6.0 * np.pi * SolarUnits.GM / (SolarUnits.c**2 * a_m * (1.0 - e**2))
        analytical_arcsec = np.degrees(analytical_rad) * 3600.0
        ax.plot(a_au, analytical_arcsec, "-", color=color, linewidth=0.8, alpha=0.7)

    ax.set_xlabel(r"$a$ (AU)")
    ax.set_ylabel(r"$\Delta\varpi$ (arcsec / orbit)")
    ax.legend(frameon=False, loc="upper right")
    _setup_axes(ax)

    if path:
        save_figure_if_changed(fig, path, dpi=150, bbox_inches="tight")
        plt.close(fig)
    else:
        plt.show()


if __name__ == "__main__":  # pragma: no cover
    plot_gr_precession()

    console = Console()
    e_mercury, a_mercury = 0.205630, 0.387098 * SolarUnits.AU
    precession_rad_mercury = compute_delta_psi(a_mercury * (1 - e_mercury), e_mercury) - 2.0 * np.pi
    precession_arcsec_mercury = np.degrees(precession_rad_mercury) * 3600.0
    orbits_per_century = 100.0 * 365.25 / 87.969
    console.print(
        f"[green]Mercury: {precession_arcsec_mercury:.4f} arcsec/orbit = "
        f"{precession_arcsec_mercury * orbits_per_century:.1f} arcsec/century[/green]"
    )
