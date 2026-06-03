"""Line-of-sight velocity dispersion of a Hernquist sphere with constant anisotropy.

Problem 6.7. The Hernquist (1990) model is self-consistent: the stars are both the
tracer and the source of the potential,

    rho(r) = M a / [2 pi r (r + a)^3],   M(<r) = M r^2 / (r + a)^2,   Phi(r) = -G M / (r + a).

For constant anisotropy beta the radial Jeans equation integrates to
``rho sigma_r^2(r) = r^(-2 beta) int_r^inf dr' r'^(2 beta) rho(r') G M(<r')/r'^2``, and the
projected dispersion follows from the Abel projection (Eq. 6.18),

    Sigma(R) sigma_los^2(R) = 2 int_R^inf dr (1 - beta R^2/r^2) rho sigma_r^2 r / sqrt(r^2 - R^2).

Both projection integrals use the substitution r = sqrt(R^2 + u^2), which removes the
1/sqrt(r^2 - R^2) singularity. Everything is computed in units G = M = a = 1, so radii are
in units of a and dispersions in units of sqrt(G M / a).
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

# Constant anisotropies to compare (tangential, isotropic, radial).
BETAS: tuple[float, ...] = (-0.5, 0.0, 0.5)
LINE_STYLES: tuple[str, ...] = ("-", "--", ":")
LINE_LABELS: tuple[str, ...] = (r"$\beta = -1/2$", r"$\beta = 0$", r"$\beta = 1/2$")


def hernquist_density(r: float) -> float:
    """Hernquist density in units G = M = a = 1."""
    return float(1.0 / (2.0 * np.pi * r * (1.0 + r) ** 3))


def hernquist_mass(r: float) -> float:
    """Hernquist enclosed mass M(<r) in units G = M = a = 1."""
    return float(r**2 / (1.0 + r) ** 2)


def surface_density(radius: float) -> float:
    """Projected surface density Sigma(R) via Abel projection (units G = M = a = 1)."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", IntegrationWarning)
        value, _ = quad(lambda u: hernquist_density(np.sqrt(radius**2 + u**2)), 0.0, np.inf, limit=200)
    return float(2.0 * value)


def _build_jeans_g(beta: float, n_grid: int) -> Callable[[float], float]:
    """Tabulate rho(r) sigma_r^2(r) on a log grid and return a log-log interpolator."""
    grid = np.logspace(-4.0, 3.0, n_grid)

    def integrand(r_prime: float) -> float:
        return float(r_prime ** (2.0 * beta) * hernquist_density(r_prime) * hernquist_mass(r_prime) / r_prime**2)

    g_values = np.empty(n_grid)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", IntegrationWarning)
        for i, r in enumerate(grid):
            value, _ = quad(integrand, r, np.inf, limit=200)
            g_values[i] = r ** (-2.0 * beta) * value

    log_g = interp1d(np.log(grid), np.log(np.clip(g_values, 1e-300, None)), fill_value="extrapolate")
    return lambda r: float(np.exp(log_g(np.log(r))))


def sigma_los_profile(radii: npt.NDArray[np.float64], beta: float, n_grid: int = 120) -> npt.NDArray[np.float64]:
    """Line-of-sight velocity dispersion profile in units of sqrt(G M / a).

    Parameters
    ----------
    radii : ndarray
        Projected radii in units of the scale radius a.
    beta : float
        Constant velocity anisotropy.
    n_grid : int, optional
        Number of radial grid points for the Jeans tabulation. Default is 120.

    Returns
    -------
    ndarray
        sigma_los / sqrt(G M / a) at each projected radius.
    """
    jeans_g = _build_jeans_g(beta, n_grid)

    def projected(radius: float) -> float:
        def integrand(u: float) -> float:
            r = np.sqrt(radius**2 + u**2)
            return float((1.0 - beta * radius**2 / r**2) * jeans_g(r))

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", IntegrationWarning)
            value, _ = quad(integrand, 0.0, np.inf, limit=200)
        return float(np.sqrt(2.0 * value / surface_density(radius)))

    return np.array([projected(float(radius)) for radius in radii])


def _setup_axes(ax: Axes) -> None:
    """Apply standard tick formatting to an axes object."""
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.tick_params(which="minor", length=3, color="gray", direction="in")
    ax.tick_params(which="major", length=6, direction="in")
    ax.tick_params(top=True, right=True, which="both")


@register_asset("p06_07_hernquist_los_dispersion.png")
def plot_hernquist_los_dispersion(path: Path | None = None) -> None:
    """Plot sigma_los/sqrt(GM/a) versus R/a for constant anisotropies.

    Parameters
    ----------
    path : Path or None
        If provided, save the figure to this path. Otherwise display it.
    """
    radii = np.logspace(np.log10(0.03), np.log10(10.0), 40)

    fig: Figure
    ax: Axes
    fig, ax = plt.subplots(figsize=(6.0, 4.6))
    fig.subplots_adjust(left=0.13, right=0.96, top=0.95, bottom=0.13)

    for beta, style, label in zip(BETAS, LINE_STYLES, LINE_LABELS):
        sigma = sigma_los_profile(radii, beta)
        ax.plot(radii, sigma, style, color="#000000", linewidth=2, label=label)

    ax.set_xscale("log")
    ax.set_xlabel(r"$R\,/\,a$")
    ax.set_ylabel(r"$\sigma_{\mathrm{los}}\,/\,\sqrt{GM/a}$")
    ax.set_xlim(0.03, 10.0)
    ax.legend(frameon=False)
    _setup_axes(ax)

    if path:
        save_figure_if_changed(fig, path, dpi=150, bbox_inches="tight")
        plt.close(fig)
    else:
        plt.show()


if __name__ == "__main__":  # pragma: no cover
    plot_hernquist_los_dispersion()
