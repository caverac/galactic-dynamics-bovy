"""Radial velocity dispersion of an NFW halo from the constant-anisotropy Jeans equation.

Problem 5.12 (part c). For a non-rotating spherical system with constant
anisotropy beta, the Jeans equation integrates to

    rho(r) sigma_r^2(r) = r^{-2 beta} int_r^inf dr' r'^{2 beta} rho(r') dPhi/dr',

with dPhi/dr' = G M(<r') / r'^2. We evaluate this integral directly (without
galpy.df.jeans) for an NFW halo

    rho(r) = rho_s / [ (r/r_s) (1 + r/r_s)^2 ],
    M(<r)  = 4 pi rho_s r_s^3 [ ln(1 + r/r_s) - (r/r_s)/(1 + r/r_s) ],

scaled to the Milky-Way-like halo with virial mass 7e11 Msun and concentration
11.5 (Delta_v = 200 times the critical density).

Units: radii in kpc, masses in Msun, dispersions in km/s, using
G = 4.300917270e-6 kpc (km/s)^2 / Msun.
"""

from dataclasses import dataclass
from pathlib import Path

from matplotlib import pyplot as plt
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.ticker import AutoMinorLocator
import numpy as np
import numpy.typing as npt
from scipy.integrate import quad

from galactic_dynamics_bovy.utils.assets import register_asset, save_figure_if_changed

# Gravitational constant in kpc (km/s)^2 / Msun.
G_KPC_KMS2_MSUN = 4.300917270e-6


@dataclass
class NFWHalo:
    """An NFW dark-matter halo specified by its scale parameters.

    Attributes
    ----------
    m_vir : float
        Virial mass in Msun.
    concentration : float
        Concentration r_vir / r_s.
    r_vir : float
        Virial radius in kpc.
    r_s : float
        Scale radius in kpc.
    rho_s : float
        Characteristic density in Msun/kpc^3.
    """

    m_vir: float
    concentration: float
    r_vir: float
    r_s: float
    rho_s: float

    @classmethod
    def from_virial_mass(
        cls,
        m_vir: float = 7e11,
        concentration: float = 11.5,
        *,
        overdensity: float = 200.0,
        hubble: float = 70.0,
    ) -> "NFWHalo":
        """Build a halo from its virial mass and concentration.

        Parameters
        ----------
        m_vir : float, optional
            Virial mass in Msun. Default is 7e11.
        concentration : float, optional
            Concentration r_vir / r_s. Default is 11.5.
        overdensity : float, optional
            Virial overdensity relative to the critical density. Default 200.
        hubble : float, optional
            Hubble constant in km/s/Mpc. Default is 70.

        Returns
        -------
        NFWHalo
            The halo with derived r_vir, r_s, and rho_s.
        """
        h0 = hubble / 1000.0  # km/s/kpc
        rho_crit = 3.0 * h0**2 / (8.0 * np.pi * G_KPC_KMS2_MSUN)
        r_vir = (3.0 * m_vir / (4.0 * np.pi * overdensity * rho_crit)) ** (1.0 / 3.0)
        r_s = r_vir / concentration
        mass_factor = np.log(1.0 + concentration) - concentration / (1.0 + concentration)
        rho_s = m_vir / (4.0 * np.pi * r_s**3 * mass_factor)
        return cls(m_vir=m_vir, concentration=concentration, r_vir=r_vir, r_s=r_s, rho_s=rho_s)

    def density(self, r: float) -> float:
        """NFW density rho(r) in Msun/kpc^3."""
        x = r / self.r_s
        return float(self.rho_s / (x * (1.0 + x) ** 2))

    def enclosed_mass(self, r: float) -> float:
        """Enclosed mass M(<r) in Msun."""
        x = r / self.r_s
        return float(4.0 * np.pi * self.rho_s * self.r_s**3 * (np.log(1.0 + x) - x / (1.0 + x)))

    def dphi_dr(self, r: float) -> float:
        """Radial gravitational acceleration dPhi/dr = G M(<r)/r^2 in (km/s)^2/kpc."""
        return G_KPC_KMS2_MSUN * self.enclosed_mass(r) / r**2

    def radial_dispersion(self, r: float, beta: float) -> float:
        """Radial velocity dispersion sigma_r(r) in km/s for constant anisotropy.

        Evaluates ``sigma_r^2 = [int_r^inf r'^{2 beta} rho dPhi/dr' dr'] /
        [rho(r) r^{2 beta}]`` by direct numerical integration.

        Parameters
        ----------
        r : float
            Radius in kpc.
        beta : float
            Constant anisotropy parameter.

        Returns
        -------
        float
            Radial velocity dispersion in km/s.
        """

        def integrand(r_prime: float) -> float:
            return float(r_prime ** (2.0 * beta) * self.density(r_prime) * self.dphi_dr(r_prime))

        integral, _ = quad(integrand, r, np.inf, limit=200)
        sigma_sq = integral / (self.density(r) * r ** (2.0 * beta))
        return float(np.sqrt(sigma_sq))


def dispersion_profile(halo: NFWHalo, radii: npt.NDArray[np.float64], beta: float) -> npt.NDArray[np.float64]:
    """Radial dispersion at each radius for a constant anisotropy beta.

    Parameters
    ----------
    halo : NFWHalo
        The halo.
    radii : ndarray
        Radii in kpc.
    beta : float
        Constant anisotropy parameter.

    Returns
    -------
    ndarray
        Radial velocity dispersions in km/s.
    """
    return np.array([halo.radial_dispersion(float(r), beta) for r in radii])


def _setup_axes(ax: Axes) -> None:
    """Apply standard tick formatting to an axes object."""
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.tick_params(which="minor", length=3, color="gray", direction="in")
    ax.tick_params(which="major", length=6, direction="in")
    ax.tick_params(top=True, right=True, which="both")


@register_asset("p05_12_nfw_dispersion.png")
def plot_nfw_dispersion(path: Path | None = None) -> None:
    """Plot the NFW radial dispersion for beta = 0 and beta = 0.5.

    Parameters
    ----------
    path : Path or None
        If provided, save the figure to this path. Otherwise display it.
    """
    halo = NFWHalo.from_virial_mass()
    radii = np.logspace(0.0, np.log10(300.0), 200)

    fig: Figure
    ax: Axes
    fig, ax = plt.subplots(figsize=(6.0, 4.6))
    fig.subplots_adjust(left=0.12, right=0.96, top=0.95, bottom=0.13)

    styles = [(0.0, "-", "#000000", r"$\beta = 0$"), (0.5, "--", "#4B4B4B", r"$\beta = 0.5$")]
    for beta, linestyle, color, label in styles:
        sigma = dispersion_profile(halo, radii, beta)
        ax.semilogx(radii, sigma, linestyle, color=color, linewidth=2, label=label)

    ax.axvline(halo.r_s, color="#737373", linestyle=":", linewidth=1.0)
    ax.text(halo.r_s, ax.get_ylim()[1], r"$r_s$", color="#737373", ha="center", va="bottom")
    ax.set_xlabel(r"$r\; (\mathrm{kpc})$")
    ax.set_ylabel(r"$\sigma_r\; (\mathrm{km\,s^{-1}})$")
    ax.set_xlim(1.0, 300.0)
    ax.legend(frameon=False)
    _setup_axes(ax)

    if path:
        save_figure_if_changed(fig, path, dpi=150, bbox_inches="tight")
        plt.close(fig)
    else:
        plt.show()


if __name__ == "__main__":  # pragma: no cover
    plot_nfw_dispersion()
