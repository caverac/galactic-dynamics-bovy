"""Local Group timing argument: the Milky Way-M31 mass from a radial Kepler orbit.

Problem 6.5. Treating the Milky Way-M31 pair as the relative (displacement)
coordinate in a point-mass potential of total mass M, the orbit is a bound
radial (eccentricity 1) Kepler orbit. Using the parametric solution from
Section 4.2.2,

    r(eta)  = a (1 - cos eta),
    t(eta)  = sqrt(a^3 / (G M)) (eta - sin eta),

the present-day radial velocity is

    v = dr/dt = sqrt(G M / a) sin(eta) / (1 - cos eta).

Eliminating a and M, the three observables (separation r, velocity v, age t)
satisfy a single transcendental equation for the eccentric anomaly eta,

    v t / r = sin(eta) (eta - sin eta) / (1 - cos eta)^2.

The pair is approaching (v < 0), so the system is past turnaround on the
infalling branch eta in (pi, 2pi). Solving for eta gives the semi-major axis
a = r / (1 - cos eta), the total mass M = a^3 (eta - sin eta)^2 / (G t^2), and
the maximum past separation (apocenter) 2a.

Units: kpc, Gyr, Msun.
"""

from dataclasses import dataclass
from pathlib import Path

from matplotlib import pyplot as plt
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.ticker import AutoMinorLocator
import numpy as np
import numpy.typing as npt
from scipy.optimize import brentq

from galactic_dynamics_bovy.utils.assets import register_asset, save_figure_if_changed
from galactic_dynamics_bovy.utils.units import GalacticUnits

# GalacticUnits expresses G in (kpc, 10^10 Msun, Gyr); divide by 1e10 for per-Msun.
# Gravitational constant in kpc^3 / (Msun Gyr^2).
G_KPC3_MSUN_GYR2 = GalacticUnits.G / 1e10
# 1 km/s in kpc/Gyr, from the ratio of G in the Gyr- and (km/s)-based unit systems:
# G [kpc^3/(M Gyr^2)] / G_kms [kpc (km/s)^2/M] = (kpc/Gyr)^2 / (km/s)^2.
KMS_TO_KPC_PER_GYR = (GalacticUnits.G / GalacticUnits.G_kms) ** 0.5


@dataclass
class TimingArgumentSolution:
    """Solution of the Local Group timing argument.

    Attributes
    ----------
    eta : float
        Present-day eccentric anomaly in radians (infalling branch).
    semi_major_axis : float
        Semi-major axis a of the radial orbit in kpc.
    mass : float
        Total Local Group mass in Msun.
    apocenter : float
        Maximum past separation 2a in kpc.
    separation : float
        Present-day separation in kpc.
    velocity : float
        Present-day relative radial velocity in km/s.
    age : float
        Age of the universe used, in Gyr.
    """

    eta: float
    semi_major_axis: float
    mass: float
    apocenter: float
    separation: float
    velocity: float
    age: float


def solve_timing_argument(
    *,
    separation: float = 740.0,
    velocity: float = -125.0,
    age: float = 13.8,
) -> TimingArgumentSolution:
    """Solve the timing argument for the Local Group mass.

    Parameters
    ----------
    separation : float, optional
        Present-day MW-M31 separation in kpc. Default is 740.
    velocity : float, optional
        Present-day relative radial velocity in km/s (negative = approaching).
        Default is -125.
    age : float, optional
        Age of the universe in Gyr. Default is 13.8.

    Returns
    -------
    TimingArgumentSolution
        The eccentric anomaly, semi-major axis, total mass, and apocenter.
    """
    v_kpc_per_gyr = velocity * KMS_TO_KPC_PER_GYR
    vtr = v_kpc_per_gyr * age / separation

    def excess(eta: float) -> float:
        return float(np.sin(eta) * (eta - np.sin(eta)) / (1.0 - np.cos(eta)) ** 2 - vtr)

    # Infalling branch: past apocenter (eta = pi), heading back toward eta = 2 pi.
    eta = float(brentq(excess, np.pi + 1e-6, 2.0 * np.pi - 1e-3))
    a = separation / (1.0 - np.cos(eta))
    mass = a**3 * (eta - np.sin(eta)) ** 2 / (G_KPC3_MSUN_GYR2 * age**2)
    return TimingArgumentSolution(
        eta=eta,
        semi_major_axis=a,
        mass=mass,
        apocenter=2.0 * a,
        separation=separation,
        velocity=velocity,
        age=age,
    )


def orbit_track(
    solution: TimingArgumentSolution, n_points: int = 400
) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]:
    """Return the separation history (time, separation) from the Big Bang to now.

    Parameters
    ----------
    solution : TimingArgumentSolution
        A solved timing argument.
    n_points : int, optional
        Number of samples along the orbit. Default is 400.

    Returns
    -------
    tuple of ndarray
        Time in Gyr (from 0 to the age) and separation in kpc.
    """
    eta = np.linspace(0.0, solution.eta, n_points)
    time_factor = np.sqrt(solution.semi_major_axis**3 / (G_KPC3_MSUN_GYR2 * solution.mass))
    time = time_factor * (eta - np.sin(eta))
    separation = solution.semi_major_axis * (1.0 - np.cos(eta))
    return time, separation


def _setup_axes(ax: Axes) -> None:
    """Apply standard tick formatting to an axes object."""
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.tick_params(which="minor", length=3, color="gray", direction="in")
    ax.tick_params(which="major", length=6, direction="in")
    ax.tick_params(top=True, right=True, which="both")


@register_asset("p06_05_local_group_timing.png")
def plot_local_group_timing(path: Path | None = None) -> None:
    """Plot the MW-M31 separation history from the timing-argument orbit.

    Parameters
    ----------
    path : Path or None
        If provided, save the figure to this path. Otherwise display it.
    """
    solution = solve_timing_argument()
    time, separation = orbit_track(solution)
    apo_time = np.sqrt(solution.semi_major_axis**3 / (G_KPC3_MSUN_GYR2 * solution.mass)) * np.pi

    fig: Figure
    ax: Axes
    fig, ax = plt.subplots(figsize=(6.0, 4.6))
    fig.subplots_adjust(left=0.13, right=0.96, top=0.95, bottom=0.13)

    ax.plot(time, separation, "-", color="#000000", linewidth=2)
    ax.plot([solution.age], [solution.separation], "o", color="#000000", markersize=7, label="today")
    ax.plot([apo_time], [solution.apocenter], "s", color="#BABABA", markersize=7, label="turnaround")
    ax.axhline(solution.apocenter, color="#737373", linestyle=":", linewidth=1.0)
    ax.set_xlabel(r"$t\; (\mathrm{Gyr})$")
    ax.set_ylabel(r"MW$-$M31 separation $r\; (\mathrm{kpc})$")
    ax.set_xlim(0.0, solution.age)
    ax.set_ylim(0.0, 1.1 * solution.apocenter)
    ax.legend(frameon=False)
    _setup_axes(ax)

    if path:
        save_figure_if_changed(fig, path, dpi=150, bbox_inches="tight")
        plt.close(fig)
    else:
        plt.show()


if __name__ == "__main__":  # pragma: no cover
    plot_local_group_timing()
