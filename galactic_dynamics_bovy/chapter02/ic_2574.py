"""IC 2574 rotation curve analysis.

This module loads and plots the rotation curve of IC 2574 from the SPARC database,
demonstrating the near-linear rise of circular velocity in the inner regions.

Data source:
    de Blok et al. (2008), "High-Resolution Rotation Curves and Galaxy Mass
    Models from THINGS", AJ, 136, 2648. DOI: 10.1088/0004-6256/136/6/2648

    SPARC Database: https://astroweb.case.edu/SPARC/
"""

from pathlib import Path

from matplotlib import pyplot as plt
from matplotlib.axes import Axes
from matplotlib.figure import Figure
import numpy as np
import numpy.typing as npt
from scipy import stats


def load_ic2574_data() -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64], npt.NDArray[np.float64]]:
    """Load IC 2574 rotation curve data from SPARC.

    Returns
    -------
    radius : ndarray
        Galactocentric radius in kpc.
    velocity : ndarray
        Observed circular velocity in km/s.
    error : ndarray
        Velocity uncertainty in km/s.
    """
    data_path = Path(__file__).parent.parent / "assets" / "IC2574_rotmod.dat"
    data = np.loadtxt(data_path, comments="#")
    radius = data[:, 0]
    velocity = data[:, 1]
    error = data[:, 2]
    return radius, velocity, error


def fit_linear_rotation(
    radius: npt.NDArray[np.float64],
    velocity: npt.NDArray[np.float64],
    r_max: float = 6.0,
) -> tuple[float, float]:
    """Fit a linear model v_c = slope * r to the inner rotation curve.

    Parameters
    ----------
    radius : ndarray
        Galactocentric radius in kpc.
    velocity : ndarray
        Observed circular velocity in km/s.
    r_max : float, optional
        Maximum radius for the fit in kpc (default: 6.0).

    Returns
    -------
    slope : float
        Best-fit slope in km/s/kpc.
    intercept : float
        Best-fit intercept in km/s.
    """
    mask = radius < r_max
    result = stats.linregress(radius[mask], velocity[mask])
    return result.slope, result.intercept


def plot_rotation_curve(path: Path | None = None) -> None:
    """Plot IC 2574 rotation curve with linear fit for r < 6 kpc.

    Parameters
    ----------
    path : Path, optional
        If provided, save the figure to this path. Otherwise, display with plt.show().
    """
    radius, velocity, error = load_ic2574_data()
    slope, intercept = fit_linear_rotation(radius, velocity, r_max=6.0)

    fig: Figure
    ax: Axes
    fig, ax = plt.subplots(1, 1, figsize=(6.5, 5))

    # Plot observed data
    ax.errorbar(
        radius,
        velocity,
        yerr=error,
        fmt="o",
        color="#000000",
        ecolor="#000000",
        alpha=1.0,
        capsize=3,
        label="Observed (THINGS)",
    )

    # Plot linear fit for r < 6 kpc
    r_fit = np.linspace(0, 6, 100)
    v_fit = slope * r_fit + intercept
    ax.plot(
        r_fit,
        v_fit,
        "--",
        color="#535353",
        linewidth=2,
        label=f"Linear fit (r < 6 kpc): $v_c = {slope:.1f} \\cdot r + {intercept:.1f}$",
    )

    ax.set_xlabel(r"$r$ (kpc)")
    ax.set_ylabel(r"$v_c$ (km/s)")
    ax.legend(loc="lower right", frameon=False)
    ax.set_xlim(0, 11)
    ax.set_ylim(0, 80)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()

    if path:
        fig.savefig(path, dpi=150, bbox_inches="tight")
        plt.close(fig)
    else:
        plt.show()
