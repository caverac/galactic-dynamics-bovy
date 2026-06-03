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
import warnings

from matplotlib import pyplot as plt
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.ticker import AutoMinorLocator
import numpy as np
import numpy.typing as npt
from rich.console import Console
from rich.table import Table
from scipy.integrate import IntegrationWarning, quad, solve_ivp
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


@dataclass
class LambdaTimingSolution:
    """Timing-argument solution including dark energy.

    Attributes
    ----------
    mass : float
        Total Local Group mass in Msun.
    apocenter : float
        Maximum past separation in kpc.
    separation : float
        Present-day separation in kpc.
    velocity : float
        Present-day relative radial velocity in km/s.
    age : float
        Age of the universe used, in Gyr.
    omega_lambda : float
        Dark-energy density parameter.
    hubble : float
        Present-day Hubble constant in km/s/Mpc.
    """

    mass: float
    apocenter: float
    separation: float
    velocity: float
    age: float
    omega_lambda: float
    hubble: float


def _lambda_acceleration_constant(omega_lambda: float, hubble: float) -> float:
    """Return the dark-energy acceleration coefficient q = Omega_Lambda H0^2 in 1/Gyr^2."""
    h0 = GalacticUnits(h=hubble / 100.0).H0
    return omega_lambda * h0**2


def solve_timing_argument_lambda(
    *,
    separation: float = 740.0,
    velocity: float = -125.0,
    age: float = 13.8,
    omega_lambda: float = 0.7,
    hubble: float = 70.0,
) -> LambdaTimingSolution:
    """Solve the timing argument with dark energy for the Local Group mass.

    The relative coordinate moves in the potential
    ``Phi(r) = -G M / r - (1/2) Omega_Lambda H0^2 r^2``. Energy conservation
    fixes M for a trial apocenter, and the age constraint (time from r=0 to
    apocenter plus apocenter to the present separation) is solved for the
    apocenter.

    Parameters
    ----------
    separation : float, optional
        Present-day separation in kpc. Default is 740.
    velocity : float, optional
        Present-day relative radial velocity in km/s. Default is -125.
    age : float, optional
        Age of the universe in Gyr. Default is 13.8.
    omega_lambda : float, optional
        Dark-energy density parameter. Default is 0.7.
    hubble : float, optional
        Present-day Hubble constant in km/s/Mpc. Default is 70.

    Returns
    -------
    LambdaTimingSolution
        The total mass and apocenter.
    """
    q = _lambda_acceleration_constant(omega_lambda, hubble)
    v_kpc_per_gyr = velocity * KMS_TO_KPC_PER_GYR

    def potential(r: float, mass: float) -> float:
        return -G_KPC3_MSUN_GYR2 * mass / r - 0.5 * q * r**2

    def mass_for_apocenter(r_apo: float) -> float:
        # Energy conservation: Phi(r_apo) = Phi(separation) + v^2/2 (v = 0 at apocenter).
        numerator = 0.5 * q * (r_apo**2 - separation**2) + 0.5 * v_kpc_per_gyr**2
        return float(numerator / ((1.0 / separation - 1.0 / r_apo) * G_KPC3_MSUN_GYR2))

    def model_age(r_apo: float) -> float:
        mass = mass_for_apocenter(r_apo)

        def integrand(r: float) -> float:
            return float(1.0 / np.sqrt(2.0 * (potential(r_apo, mass) - potential(r, mass))))

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", IntegrationWarning)
            time_to_apo, _ = quad(integrand, 0.0, r_apo, limit=200)
            time_from_apo, _ = quad(integrand, separation, r_apo, limit=200)
        return time_to_apo + time_from_apo

    r_apo = float(brentq(lambda r: model_age(r) - age, separation + 1.0, 1600.0))
    return LambdaTimingSolution(
        mass=mass_for_apocenter(r_apo),
        apocenter=r_apo,
        separation=separation,
        velocity=velocity,
        age=age,
        omega_lambda=omega_lambda,
        hubble=hubble,
    )


def lambda_orbit_track(
    solution: LambdaTimingSolution, n_points: int = 400
) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]:
    """Integrate the dark-energy orbit to get the separation history.

    Parameters
    ----------
    solution : LambdaTimingSolution
        A solved dark-energy timing argument.
    n_points : int, optional
        Number of samples. Default is 400.

    Returns
    -------
    tuple of ndarray
        Time in Gyr (from 0 to the age) and separation in kpc.
    """
    q = _lambda_acceleration_constant(solution.omega_lambda, solution.hubble)
    mass = solution.mass
    v_kpc_per_gyr = solution.velocity * KMS_TO_KPC_PER_GYR
    energy = -G_KPC3_MSUN_GYR2 * mass / solution.separation - 0.5 * q * solution.separation**2 + 0.5 * v_kpc_per_gyr**2

    r_start = 0.5
    v_start = np.sqrt(2.0 * (energy + G_KPC3_MSUN_GYR2 * mass / r_start + 0.5 * q * r_start**2))

    def rhs(_t: float, state: npt.NDArray[np.float64]) -> list[float]:
        r, v = state
        return [v, -G_KPC3_MSUN_GYR2 * mass / r**2 + q * r]

    sol = solve_ivp(
        rhs,
        (0.0, solution.age),
        np.array([r_start, v_start]),
        t_eval=np.linspace(0.0, solution.age, n_points),
        rtol=1e-10,
        atol=1e-10,
    )
    return sol.t, sol.y[0]


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
    matter = solve_timing_argument()
    time_m, separation_m = orbit_track(matter)
    lam = solve_timing_argument_lambda()
    time_l, separation_l = lambda_orbit_track(lam)

    fig: Figure
    ax: Axes
    fig, ax = plt.subplots(figsize=(6.0, 4.6))
    fig.subplots_adjust(left=0.13, right=0.96, top=0.95, bottom=0.13)

    ax.plot(time_m, separation_m, "-", color="#000000", linewidth=2, label="matter only")
    ax.plot(time_l, separation_l, "--", color="#AA3377", linewidth=2, label=r"$\Lambda$CDM ($\Omega_\Lambda=0.7$)")
    # Mark each apocenter with a dot (not added to the legend).
    ax.plot(time_m[np.argmax(separation_m)], np.max(separation_m), "o", color="#000000", markersize=6)
    ax.plot(time_l[np.argmax(separation_l)], np.max(separation_l), "o", color="#AA3377", markersize=6)
    ax.set_xlabel(r"$t\; (\mathrm{Gyr})$")
    ax.set_ylabel(r"MW$-$M31 separation $r\; (\mathrm{kpc})$")
    ax.set_xlim(0.0, matter.age)
    ax.set_ylim(0.0, 1.1 * max(matter.apocenter, lam.apocenter))
    ax.legend(frameon=False)
    _setup_axes(ax)

    if path:
        save_figure_if_changed(fig, path, dpi=150, bbox_inches="tight")
        plt.close(fig)
    else:
        plt.show()


if __name__ == "__main__":  # pragma: no cover
    console = Console()
    matter_sol = solve_timing_argument()
    lambda_sol = solve_timing_argument_lambda()
    vt_ratio = matter_sol.velocity * KMS_TO_KPC_PER_GYR * matter_sol.age / matter_sol.separation

    table = Table(title="Problem 6.5: Local Group timing argument (r = 740 kpc, v = -125 km/s, age = 13.8 Gyr)")
    table.add_column("Quantity")
    table.add_column("Matter only", justify="right")
    table.add_column(f"Lambda CDM (Omega_L = {lambda_sol.omega_lambda})", justify="right")
    table.add_row("v t / r", f"{vt_ratio:.3f}", "-")
    table.add_row("eta (rad)", f"{matter_sol.eta:.3f}", "-")
    table.add_row("semi-major axis a (kpc)", f"{matter_sol.semi_major_axis:.1f}", "-")
    table.add_row("mass M (Msun)", f"{matter_sol.mass:.3e}", f"{lambda_sol.mass:.3e}")
    table.add_row("apocenter (kpc)", f"{matter_sol.apocenter:.1f}", f"{lambda_sol.apocenter:.1f}")
    console.print(table)
    plot_local_group_timing()
