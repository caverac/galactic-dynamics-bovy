"""Circular velocity profiles for Jaffe, Hernquist, and NFW models.

This module computes and plots the dimensionless circular velocity
v_c / sqrt(4 pi G rho_0 a^2) as a function of r/a for three common
spherical density profiles.

The three profiles are:

- Jaffe: rho = rho_0 / [(r/a)^2 (1 + r/a)^2]
- Hernquist: rho = rho_0 / [(r/a) (1 + r/a)^3]
- NFW: rho = rho_0 / [(r/a) (1 + r/a)^2]
"""

from pathlib import Path

from matplotlib import pyplot as plt
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.ticker import AutoMinorLocator
import numpy as np
import numpy.typing as npt


def vcirc_jaffe(x: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
    """Compute normalized circular velocity for the Jaffe profile.

    Parameters
    ----------
    x : ndarray
        Dimensionless radius r/a.

    Returns
    -------
    ndarray
        Normalized circular velocity v_c / sqrt(4 pi G rho_0 a^2).
    """
    return 1.0 / np.sqrt(1.0 + x)


def vcirc_hernquist(x: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
    """Compute normalized circular velocity for the Hernquist profile.

    Parameters
    ----------
    x : ndarray
        Dimensionless radius r/a.

    Returns
    -------
    ndarray
        Normalized circular velocity v_c / sqrt(4 pi G rho_0 a^2).
    """
    return np.sqrt(x / (2.0 * (1.0 + x) ** 2))


def vcirc_nfw(x: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
    """Compute normalized circular velocity for the NFW profile.

    Parameters
    ----------
    x : ndarray
        Dimensionless radius r/a.

    Returns
    -------
    ndarray
        Normalized circular velocity v_c / sqrt(4 pi G rho_0 a^2).
    """
    f_x = np.log(1.0 + x) - x / (1.0 + x)
    return np.sqrt(f_x / x)


def plot_vcirc_profiles(path: Path | None = None) -> None:
    """Plot circular velocity curves for Jaffe, Hernquist, and NFW profiles.

    Creates a plot showing the dimensionless circular velocity
    v_c / sqrt(4 pi G rho_0 a^2) as a function of r/a for three profiles.

    Parameters
    ----------
    path : Path or None
        If provided, save the figure to this path. Otherwise display it.
    """
    x = np.linspace(1e-3, 6, 500)

    fig: Figure
    ax: Axes
    fig, ax = plt.subplots(figsize=(6.5, 4.5))

    ax.plot(x, vcirc_jaffe(x), "-", color="#000000", linewidth=2, label="Jaffe")
    ax.plot(x, vcirc_hernquist(x), "--", color="#303030", linewidth=2, label="Hernquist")
    ax.plot(x, vcirc_nfw(x), "-.", color="#737373", linewidth=2, label="NFW")

    ax.set_xlabel(r"$r/a$")
    ax.set_ylabel(r"$v_c / \sqrt{4\pi G \rho_0 a^2}$")
    ax.set_xlim(0.0, 6)
    ax.set_ylim(0, 1.1)

    # Only set minor locator on y-axis (x-axis is log scale)
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.tick_params(which="minor", length=3, color="gray", direction="in")
    ax.tick_params(which="major", length=6, direction="in")
    ax.tick_params(top=True, right=True, which="both")

    ax.legend(loc="upper right", frameon=False)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()

    if path:
        fig.savefig(path, dpi=150, bbox_inches="tight")
        plt.close(fig)
    else:
        plt.show()
