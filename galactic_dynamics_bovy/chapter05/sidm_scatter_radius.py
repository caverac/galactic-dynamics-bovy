"""Self-interacting dark matter (SIDM) once-scattered radius for an NFW halo.

Problem 5.12 (part d). In SIDM, dark-matter particles thermalize in the dense
inner halo, where the profile becomes the cored isothermal sphere of parts a-b.
The boundary with the outer (collisionless) NFW profile is the radius r_1 where
a particle has scattered on average once over the halo age. The scattering rate
per particle is

    Gamma(r) = (sigma/m) (4/sqrt(pi)) sigma_r(r) rho(r),

where sigma/m is the SIDM cross section per unit mass, sigma_r(r) the radial
velocity dispersion (Problem 5.12c), and rho(r) the NFW density. The radius r_1
solves Gamma(r_1) * t_age = 1.

We use the NFW halo and constant-anisotropy dispersions from
``nfw_jeans_dispersion`` and report r_1 for beta = 0 and beta = 0.5 with
sigma/m = 1 cm^2/g and t_age = 10 Gyr.
"""

from pathlib import Path

from matplotlib import pyplot as plt
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.ticker import AutoMinorLocator
import numpy as np
import numpy.typing as npt
from scipy.optimize import brentq

from galactic_dynamics_bovy.chapter05.nfw_jeans_dispersion import NFWHalo
from galactic_dynamics_bovy.utils.assets import register_asset, save_figure_if_changed

# Unit conversions to make Gamma * t_age dimensionless when sigma/m is in cm^2/g,
# sigma_r in km/s, rho in Msun/kpc^3, and t_age in Gyr.
_MSUN_IN_G = 1.988409e33
_KPC_IN_CM = 3.0856775814913673e21
_GYR_IN_S = 3.15576e16
_KM_S_IN_CM_S = 1.0e5

# Gamma * t_age = SCATTER_PREFACTOR * (sigma/m) * t_age[Gyr] * sigma_r[km/s] * rho[Msun/kpc^3].
SCATTER_PREFACTOR = (4.0 / np.sqrt(np.pi)) * _KM_S_IN_CM_S * (_MSUN_IN_G / _KPC_IN_CM**3) * _GYR_IN_S

# Cross-section values for the bonus plot, in cm^2/g.
CROSS_SECTIONS: tuple[float, ...] = (0.5, 1.0, 5.0, 10.0)


def expected_scatters(
    halo: NFWHalo,
    r: float,
    beta: float,
    *,
    cross_section: float = 1.0,
    age: float = 10.0,
) -> float:
    """Compute the expected number of scatters Gamma(r) * t_age at radius r.

    Parameters
    ----------
    halo : NFWHalo
        The dark-matter halo.
    r : float
        Radius in kpc.
    beta : float
        Constant anisotropy parameter.
    cross_section : float, optional
        SIDM cross section per unit mass sigma/m in cm^2/g. Default is 1.0.
    age : float, optional
        Halo age in Gyr. Default is 10.

    Returns
    -------
    float
        Dimensionless expected number of scatters.
    """
    sigma_r = halo.radial_dispersion(r, beta)
    return float(SCATTER_PREFACTOR * cross_section * age * sigma_r * halo.density(r))


def scatter_radius(
    halo: NFWHalo,
    beta: float,
    *,
    cross_section: float = 1.0,
    age: float = 10.0,
    r_min: float = 0.3,
    r_max: float = 300.0,
) -> float:
    """Radius r_1 (kpc) where a particle has scattered once over the halo age.

    Solves ``expected_scatters(r_1) = 1`` by bracketed root finding.

    Parameters
    ----------
    halo : NFWHalo
        The dark-matter halo.
    beta : float
        Constant anisotropy parameter.
    cross_section : float, optional
        SIDM cross section per unit mass in cm^2/g. Default is 1.0.
    age : float, optional
        Halo age in Gyr. Default is 10.
    r_min, r_max : float, optional
        Bracket for the root in kpc. Default 0.3 to 300.

    Returns
    -------
    float
        The once-scattered radius r_1 in kpc.
    """

    def excess(r: float) -> float:
        return expected_scatters(halo, r, beta, cross_section=cross_section, age=age) - 1.0

    return float(brentq(excess, r_min, r_max))


def scatter_radius_vs_beta(
    halo: NFWHalo,
    betas: npt.NDArray[np.float64],
    cross_sections: tuple[float, ...],
    *,
    age: float = 10.0,
    n_grid: int = 60,
    r_min: float = 0.3,
    r_max: float = 300.0,
) -> npt.NDArray[np.float64]:
    """Once-scattered radius r_1 over a grid of anisotropies and cross sections.

    For each beta the product ``sigma_r(r) rho(r)`` is tabulated once on a
    logarithmic radius grid and reused for every cross section (the dispersion
    integral does not depend on sigma/m), which keeps the computation fast.

    Parameters
    ----------
    halo : NFWHalo
        The dark-matter halo.
    betas : ndarray
        Anisotropy values.
    cross_sections : tuple of float
        Cross sections per unit mass in cm^2/g.
    age : float, optional
        Halo age in Gyr. Default is 10.
    n_grid : int, optional
        Number of radius grid points per beta. Default is 60.
    r_min, r_max : float, optional
        Radius grid range in kpc. Default 0.3 to 300.

    Returns
    -------
    ndarray
        Array of shape ``(len(cross_sections), len(betas))`` of r_1 in kpc.
    """
    r_grid = np.logspace(np.log10(r_min), np.log10(r_max), n_grid)
    log_r = np.log(r_grid)
    result = np.empty((len(cross_sections), len(betas)))
    for j, beta in enumerate(betas):
        product = np.array([halo.radial_dispersion(float(r), float(beta)) * halo.density(float(r)) for r in r_grid])
        # sigma_r * rho decreases outward; reverse so the interpolation abscissa increases.
        log_product = np.log(product)[::-1]
        log_r_rev = log_r[::-1]
        for i, cross_section in enumerate(cross_sections):
            target = np.log(1.0 / (SCATTER_PREFACTOR * cross_section * age))
            result[i, j] = np.exp(np.interp(target, log_product, log_r_rev))
    return result


def _setup_axes(ax: Axes) -> None:
    """Apply standard tick formatting to an axes object."""
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.tick_params(which="minor", length=3, color="gray", direction="in")
    ax.tick_params(which="major", length=6, direction="in")
    ax.tick_params(top=True, right=True, which="both")


@register_asset("p05_12_sidm_scatter_radius.png")
def plot_scatter_radius_vs_beta(path: Path | None = None) -> None:
    """Plot the once-scattered radius r_1 versus beta for several cross sections.

    Parameters
    ----------
    path : Path or None
        If provided, save the figure to this path. Otherwise display it.
    """
    halo = NFWHalo.from_virial_mass()
    betas = np.linspace(-0.5, 0.9, 30)
    radii = scatter_radius_vs_beta(halo, betas, CROSS_SECTIONS)

    fig: Figure
    ax: Axes
    fig, ax = plt.subplots(figsize=(6.0, 4.6))
    fig.subplots_adjust(left=0.12, right=0.96, top=0.95, bottom=0.13)

    grays = ["#000000", "#555555", "#888888", "#BBBBBB"]
    for cross_section, color in zip(CROSS_SECTIONS, grays):
        ax.plot(
            betas,
            radii[CROSS_SECTIONS.index(cross_section)],
            "-",
            color=color,
            linewidth=2,
            label=rf"$\sigma/m = {cross_section:g}\;\mathrm{{cm^2/g}}$",
        )

    ax.set_xlabel(r"$\beta$")
    ax.set_ylabel(r"$r_1\; (\mathrm{kpc})$")
    ax.legend(frameon=False, fontsize=9)
    _setup_axes(ax)

    if path:
        save_figure_if_changed(fig, path, dpi=150, bbox_inches="tight")
        plt.close(fig)
    else:
        plt.show()


if __name__ == "__main__":  # pragma: no cover
    plot_scatter_radius_vs_beta()
