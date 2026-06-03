"""Mass-anisotropy and cusp-core degeneracy in dwarf spheroidal Jeans fits.

Problem 6.6 (cf. Figure 6.4). Line-of-sight velocity dispersion profiles of dwarf
spheroidals are fit with the spherical Jeans equation. We investigate how the
assumed velocity anisotropy beta affects the inferred inner density slope, and
whether a cored profile can fit the data as well as a cusp.

Setup
-----
Tracer stars follow a Plummer profile with scale ``b``,

    nu(r) = (1 + r^2/b^2)^(-5/2),    Sigma(R) = (1 + R^2/b^2)^(-2),

(normalization cancels in sigma_los). The dark-matter halo is a generalized
NFW/Zhao profile with inner slope gamma and outer slope 3,

    rho(r) = rho_s (r/r_s)^(-gamma) (1 + r/r_s)^(gamma - 3).

For constant anisotropy beta the Jeans equation gives
``nu sigma_r^2(r) = r^(-2 beta) int_r^inf dr' r'^(2 beta) nu(r') G M(<r')/r'^2``,
and the projected dispersion follows from the Abel-type projection (Eq. 6.18),

    Sigma(R) sigma_los^2(R) = 2 int_R^inf dr (1 - beta R^2/r^2) nu sigma_r^2 r / sqrt(r^2 - R^2).

The projection singularity at r = R is removed with the substitution
r = sqrt(R^2 + u^2), giving ``r dr / sqrt(r^2 - R^2) = du``.

Units: kpc, km/s, Msun.
"""

from collections.abc import Callable
from pathlib import Path
import warnings

from matplotlib import pyplot as plt
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.ticker import AutoMinorLocator
import numpy as np
import numpy.typing as npt
from scipy.integrate import IntegrationWarning, quad
from scipy.interpolate import interp1d

from galactic_dynamics_bovy.utils.assets import register_asset, save_figure_if_changed
from galactic_dynamics_bovy.utils.units import GalacticUnits

# Gravitational constant in kpc (km/s)^2 / Msun.
G_KPC_KMS2_MSUN = GalacticUnits.G_kms / 1e10

# Fiducial dwarf-spheroidal scales.
TRACER_SCALE = 0.25  # Plummer scale b in kpc
HALO_SCALE = 1.0  # halo scale r_s in kpc
R_NORM = 0.5  # radius at which the models share the same enclosed mass, kpc
FIDUCIAL_RHO_S = 6.0e7  # characteristic density of the fiducial cusp, Msun/kpc^3


def _plummer_nu(r: float, tracer_scale: float) -> float:
    """Unnormalized Plummer 3D tracer density."""
    return float((1.0 + r**2 / tracer_scale**2) ** -2.5)


def _plummer_sigma(radius: float, tracer_scale: float) -> float:
    """Unnormalized Plummer surface density."""
    return float((1.0 + radius**2 / tracer_scale**2) ** -2.0)


def enclosed_mass(r: float, *, rho_s: float, gamma: float, halo_scale: float = HALO_SCALE) -> float:
    """Halo mass within radius r for a generalized NFW profile.

    Parameters
    ----------
    r : float
        Radius in kpc.
    rho_s : float
        Characteristic density in Msun/kpc^3.
    gamma : float
        Inner density slope (1 = NFW cusp, 0 = core).
    halo_scale : float, optional
        Scale radius r_s in kpc. Default is HALO_SCALE.

    Returns
    -------
    float
        Enclosed mass in Msun.
    """

    def integrand(x: float) -> float:
        return float((x / halo_scale) ** (-gamma) * (1.0 + x / halo_scale) ** (gamma - 3.0) * x**2)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", IntegrationWarning)
        value, _ = quad(integrand, 0.0, r, limit=100)
    return float(4.0 * np.pi * rho_s * value)


def density_normalization(*, target_mass: float, r_norm: float, gamma: float, halo_scale: float = HALO_SCALE) -> float:
    """Return rho_s so that the halo encloses target_mass within r_norm."""
    unit_mass = enclosed_mass(r_norm, rho_s=1.0, gamma=gamma, halo_scale=halo_scale)
    return target_mass / unit_mass


def _build_jeans_g(
    rho_s: float,
    *,
    gamma: float,
    beta: float,
    tracer_scale: float,
    halo_scale: float,
    n_grid: int,
) -> Callable[[float], float]:
    """Tabulate nu(r) sigma_r^2(r) on a log grid and return a log-log interpolator."""
    grid = np.logspace(-3.0, 1.7, n_grid)

    masses = np.array([enclosed_mass(r, rho_s=rho_s, gamma=gamma, halo_scale=halo_scale) for r in grid])
    log_mass = interp1d(np.log(grid), np.log(masses), fill_value="extrapolate")

    def mass_of(r: float) -> float:
        return float(np.exp(log_mass(np.log(r))))

    def jeans_integrand(r_prime: float) -> float:
        weight = r_prime ** (2.0 * beta) * _plummer_nu(r_prime, tracer_scale)
        return float(weight * G_KPC_KMS2_MSUN * mass_of(r_prime) / r_prime**2)

    g_values = np.empty(n_grid)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", IntegrationWarning)
        for i, r in enumerate(grid):
            value, _ = quad(jeans_integrand, r, 50.0, limit=100)
            g_values[i] = r ** (-2.0 * beta) * value

    log_g = interp1d(np.log(grid), np.log(np.clip(g_values, 1e-300, None)), fill_value="extrapolate")
    return lambda r: float(np.exp(log_g(np.log(r))))


def sigma_los_profile(
    radii: npt.NDArray[np.float64],
    *,
    rho_s: float,
    gamma: float,
    beta: float,
    tracer_scale: float = TRACER_SCALE,
    halo_scale: float = HALO_SCALE,
    n_grid: int = 80,
) -> npt.NDArray[np.float64]:
    """Projected line-of-sight velocity dispersion profile sigma_los(R).

    Parameters
    ----------
    radii : ndarray
        Projected radii in kpc.
    rho_s : float
        Halo characteristic density in Msun/kpc^3.
    gamma : float
        Inner density slope.
    beta : float
        Constant velocity anisotropy.
    tracer_scale : float, optional
        Plummer scale b in kpc. Default is TRACER_SCALE.
    halo_scale : float, optional
        Halo scale r_s in kpc. Default is HALO_SCALE.
    n_grid : int, optional
        Number of radial grid points for the Jeans tabulation. Default is 80.

    Returns
    -------
    ndarray
        Line-of-sight velocity dispersion in km/s at each projected radius.
    """
    jeans_g = _build_jeans_g(
        rho_s, gamma=gamma, beta=beta, tracer_scale=tracer_scale, halo_scale=halo_scale, n_grid=n_grid
    )

    def projected(radius: float) -> float:
        def integrand(u: float) -> float:
            r = np.sqrt(radius**2 + u**2)
            return float((1.0 - beta * radius**2 / r**2) * jeans_g(r))

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", IntegrationWarning)
            value, _ = quad(integrand, 0.0, 8.0, limit=200)
        return float(np.sqrt(2.0 * value / _plummer_sigma(radius, tracer_scale)))

    return np.array([projected(float(radius)) for radius in radii])


def _setup_axes(ax: Axes) -> None:
    """Apply standard tick formatting to an axes object."""
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.tick_params(which="minor", length=3, color="gray", direction="in")
    ax.tick_params(which="major", length=6, direction="in")
    ax.tick_params(top=True, right=True, which="both")


@register_asset("p06_06_dsph_mass_anisotropy.png")
def plot_dsph_mass_anisotropy(path: Path | None = None) -> None:
    """Plot the cusp-core / mass-anisotropy degeneracy in sigma_los(R).

    Parameters
    ----------
    path : Path or None
        If provided, save the figure to this path. Otherwise display it.
    """
    radii = np.logspace(np.log10(0.03), np.log10(10), 25)
    target_mass = enclosed_mass(R_NORM, rho_s=FIDUCIAL_RHO_S, gamma=1.0)
    rho_core = density_normalization(target_mass=target_mass, r_norm=R_NORM, gamma=0.0)

    cusp = sigma_los_profile(radii, rho_s=FIDUCIAL_RHO_S, gamma=1.0, beta=0.0)
    core_iso = sigma_los_profile(radii, rho_s=rho_core, gamma=0.0, beta=0.0)
    core_rad = sigma_los_profile(radii, rho_s=rho_core, gamma=0.0, beta=0.5)

    fig: Figure
    ax: Axes
    fig, ax = plt.subplots(figsize=(6.0, 4.6))
    fig.subplots_adjust(left=0.13, right=0.96, top=0.95, bottom=0.13)

    ax.plot(radii, cusp, "-", color="#000000", linewidth=2, label=r"cusp $\gamma=1$, $\beta=0$")
    ax.plot(radii, core_iso, "--", color="#000000", linewidth=2, label=r"core $\gamma=0$, $\beta=0$")
    ax.plot(radii, core_rad, ":", color="#000000", linewidth=2.4, label=r"core $\gamma=0$, $\beta=0.5$")
    ax.set_xscale("log")
    ax.set_xlabel(r"$R\; (\mathrm{kpc})$")
    ax.set_ylabel(r"$\sigma_{\mathrm{los}}\; (\mathrm{km\,s^{-1}})$")
    ax.set_xlim(0.03, 10)
    ax.legend(frameon=False)
    _setup_axes(ax)

    if path:
        save_figure_if_changed(fig, path, dpi=150, bbox_inches="tight")
        plt.close(fig)
    else:
        plt.show()


if __name__ == "__main__":  # pragma: no cover
    plot_dsph_mass_anisotropy()
