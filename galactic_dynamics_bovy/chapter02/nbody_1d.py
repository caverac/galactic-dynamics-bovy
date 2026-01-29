"""One-dimensional N-body gravity simulation.

This module implements an O(N log N) algorithm for computing gravitational
forces in a 1D system of N point masses.

In 1D, the force from a sheet of mass A at position x' on a test particle
at position x is:

    F = -2 pi G A sign(x - x')

For N particles of equal mass m = A/N, the force on particle i is:

    F_i = 2 pi G m (N_i^+ - N_i^-)

where N_i^+ and N_i^- are the number of particles to the right and left
of particle i, respectively.
"""

from pathlib import Path

from matplotlib import pyplot as plt
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.ticker import AutoMinorLocator
import numpy as np
import numpy.typing as npt


def compute_forces_1d(
    x: npt.NDArray[np.float64],
    m: float = 1.0,
    G: float = 1.0,
) -> npt.NDArray[np.float64]:
    """Compute gravitational forces on N particles in 1D.

    Uses an O(N log N) algorithm based on sorting.

    Parameters
    ----------
    x : ndarray
        Positions of N particles.
    m : float
        Mass of each particle (assumed equal).
    G : float
        Gravitational constant.

    Returns
    -------
    ndarray
        Force on each particle (same order as input positions).
    """
    n = len(x)
    # Get sorted indices
    sorted_indices = np.argsort(x)

    # For particle at sorted position k: N_left = k, N_right = n - 1 - k
    # F = 2*pi*G*m * (N_right - N_left) = 2*pi*G*m * (n - 1 - 2*k)
    k = np.arange(n)
    forces_sorted = 2 * np.pi * G * m * (n - 1 - 2 * k)

    # Map forces back to original order
    forces = np.empty(n)
    forces[sorted_indices] = forces_sorted

    return forces


def force_theory_uniform(
    x: npt.NDArray[np.float64],
    total_mass: float = 1.0,
    G: float = 1.0,
) -> npt.NDArray[np.float64]:
    """Compute theoretical force for uniform distribution on [-1/2, 1/2].

    For a uniform 1D mass distribution with total mass A on [-1/2, 1/2],
    the force at position x is:

        F(x) = -4 pi G A x

    Parameters
    ----------
    x : ndarray
        Positions where to evaluate the force.
    total_mass : float
        Total mass of the distribution.
    G : float
        Gravitational constant.

    Returns
    -------
    ndarray
        Theoretical force at each position.
    """
    return -4 * np.pi * G * total_mass * x


def plot_force_comparison(
    n_particles: int = 10000,
    path: Path | None = None,
) -> None:
    """Plot comparison of numerical vs theoretical 1D gravitational force.

    Creates a plot showing ln(F_numerical / F_theory) vs position for
    particles uniformly distributed on [-1/2, 1/2].

    Parameters
    ----------
    n_particles : int
        Number of particles in the simulation.
    path : Path or None
        If provided, save the figure to this path. Otherwise display it.
    """
    # Generate uniform distribution
    rng = np.random.default_rng(42)
    x = rng.uniform(-0.5, 0.5, n_particles)

    # Compute forces
    total_mass = n_particles  # Each particle has mass 1
    f_numerical = compute_forces_1d(x, m=1.0, G=1.0)
    f_theory = force_theory_uniform(x, total_mass=total_mass, G=1.0)

    # Compute error
    x_plot = x
    error = (f_numerical - f_theory) / (4 * np.pi * total_mass)

    # Sort for cleaner plotting
    sort_idx = np.argsort(x_plot)
    x_plot = x_plot[sort_idx]
    error = error[sort_idx]

    fig: Figure
    ax: Axes
    fig, ax = plt.subplots(figsize=(6.5, 4.5))

    ax.plot(x_plot, error, ".", color="#000000", markersize=1, alpha=0.5)
    ax.axhline(0, color="#888888", linestyle="--", linewidth=1)

    ax.set_xlabel(r"$x_i$")
    ax.set_ylabel(r"$(F_{\mathrm{numerical}} - F_{\mathrm{theory}}) / (4 \pi G M)$")
    ax.set_xlim(-0.5, 0.5)

    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.tick_params(which="minor", length=3, color="gray", direction="in")
    ax.tick_params(which="major", length=6, direction="in")
    ax.tick_params(top=True, right=True, which="both")

    ax.grid(True, alpha=0.3)

    plt.tight_layout()

    if path:
        fig.savefig(path, dpi=150, bbox_inches="tight")
        plt.close(fig)
    else:
        plt.show()
