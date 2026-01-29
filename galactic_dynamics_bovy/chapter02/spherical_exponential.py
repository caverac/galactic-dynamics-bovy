"""Circular velocity for spherical exponential disk and true exponential disk.

This module computes and compares the circular velocity curves for:

1. Spherical exponential disk: A spherical distribution with the same
   enclosed mass profile as an exponential disk.

2. True exponential disk: The actual rotation curve of a razor-thin
   exponential disk, involving modified Bessel functions.

The spherical exponential disk has:
    v_c^2(r) = (GM/r) [1 - (1 + r/R_d) exp(-r/R_d)]

The true exponential disk has (Freeman 1970):
    v_c^2(R) = 4 pi G Sigma_0 R_d y^2 [I_0(y)K_0(y) - I_1(y)K_1(y)]
where y = R/(2 R_d).
"""

from pathlib import Path

from matplotlib import pyplot as plt
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.ticker import AutoMinorLocator
import numpy as np
import numpy.typing as npt
from scipy.special import i0, i1, k0, k1  # pylint: disable=no-name-in-module


def vcirc_spherical_exponential(x: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
    """Compute normalized circular velocity for spherical exponential disk.

    Parameters
    ----------
    x : ndarray
        Dimensionless radius r/R_d.

    Returns
    -------
    ndarray
        Normalized circular velocity squared v_c^2 / (GM/R_d).
    """
    return (1.0 / x) * (1.0 - (1.0 + x) * np.exp(-x))


def vcirc_exponential_disk(x: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
    """Compute normalized circular velocity for true exponential disk.

    Uses the Freeman (1970) formula involving modified Bessel functions.

    Parameters
    ----------
    x : ndarray
        Dimensionless radius R/R_d.

    Returns
    -------
    ndarray
        Normalized circular velocity squared v_c^2 / (GM/R_d).
    """
    y = x / 2.0
    return 2.0 * y**2 * (i0(y) * k0(y) - i1(y) * k1(y))


def plot_exponential_vcirc(path: Path | None = None) -> None:
    """Plot circular velocity for spherical and true exponential disks.

    Creates a plot comparing the dimensionless circular velocity
    v_c / sqrt(GM/R_d) as a function of R/R_d for both profiles.

    Parameters
    ----------
    path : Path or None
        If provided, save the figure to this path. Otherwise display it.
    """
    x = np.linspace(0.01, 8, 500)

    vc2_spherical = vcirc_spherical_exponential(x)
    vc2_disk = vcirc_exponential_disk(x)

    fig: Figure
    ax: Axes
    fig, ax = plt.subplots(figsize=(6.5, 4.5))

    ax.plot(
        x,
        np.sqrt(vc2_spherical),
        "-",
        color="#000000",
        linewidth=2,
        label="Spherical exponential",
    )
    ax.plot(
        x,
        np.sqrt(vc2_disk),
        "--",
        color="#505050",
        linewidth=2,
        label="True exponential disk",
    )

    ax.set_xlabel(r"$R/R_d$")
    ax.set_ylabel(r"$v_c / \sqrt{GM/R_d}$")
    ax.set_xlim(0, 8)
    ax.set_ylim(0, 0.8)

    ax.xaxis.set_minor_locator(AutoMinorLocator())
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
