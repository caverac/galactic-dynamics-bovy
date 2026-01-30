"""Effective radius as a function of gamma for spherical systems.

This module computes the effective (half-light) radius for a family of
spherical density profiles parameterized by gamma, using numerical
integration and interpolation.

The key integrals are:

    s(x) = (3 - gamma) integral_x^inf dy y^(1-gamma) / ((1 + y)^(4-gamma) sqrt(y^2 - x^2))

    m(x) = 2*pi integral_0^x dx' s(x') x'

The effective radius x_e satisfies m(x_e) = 1/2.
"""

from pathlib import Path

from matplotlib import pyplot as plt
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.ticker import AutoMinorLocator
import numpy as np
import numpy.typing as npt
from scipy import integrate, interpolate

from galactic_dynamics_bovy.utils.assets import register_asset, save_figure_if_changed


def compute_s_direct(x: float, gamma: float) -> float:
    """Compute s(x) by direct numerical integration.

    Parameters
    ----------
    x : float
        Position at which to evaluate s.
    gamma : float
        Profile parameter.

    Returns
    -------
    float
        Value of s(x).
    """
    if x <= 0:
        # Handle x=0 case: integral from 0 to inf
        result, _ = integrate.quad(
            lambda y: y ** (1 - gamma) / (1 + y) ** (4 - gamma) / y,
            0,
            np.inf,
            limit=200,
        )
        return (3 - gamma) * result

    # For x > 0, use substitution y = x * cosh(t) to handle the singularity
    # dy = x * sinh(t) dt, sqrt(y^2 - x^2) = x * sinh(t)
    # So the sqrt cancels with dy, leaving:
    # integral = int_0^inf dt * (x*cosh(t))^(1-gamma) / (1 + x*cosh(t))^(4-gamma)

    def transformed_integrand(t: float) -> float:
        if t > 700:  # Prevent cosh overflow
            return 0.0
        y = x * np.cosh(t)
        # Use log-space to avoid overflow
        log_num = (1 - gamma) * np.log(y)
        log_den = (4 - gamma) * np.log(1 + y)
        log_result = log_num - log_den
        if log_result < -700:
            return 0.0
        return float(np.exp(log_result))

    result, _ = integrate.quad(transformed_integrand, 0, np.inf, limit=200)
    return (3 - gamma) * result


class SInterpolator:
    """Interpolator for s(x) = (3-gamma) integral_x^inf dy y^(1-gamma) / ((1+y)^(4-gamma) sqrt(y^2-x^2)).

    Parameters
    ----------
    gamma : float
        Profile parameter.
    x_min : float
        Minimum x for interpolation grid.
    x_max : float
        Maximum x for interpolation grid.
    n_points : int
        Number of grid points.
    """

    def __init__(
        self,
        gamma: float,
        x_min: float | None = None,
        x_max: float = 1e4,
        n_points: int = 200,
    ) -> None:
        """Initialize the interpolator for s(x)."""
        self.gamma = gamma

        if x_min is None:
            x_min = 1e-9
        self.x_min = x_min
        self.x_max = x_max

        # Build logarithmic grid
        self._log_x = np.linspace(np.log10(x_min), np.log10(x_max), n_points)
        self._x_grid = 10**self._log_x

        # Compute s on the grid
        self._s_grid = np.array([compute_s_direct(x, gamma) for x in self._x_grid])

        # Build interpolator (in log-log space for better accuracy)
        self._log_s = np.log10(np.maximum(self._s_grid, 1e-100))
        self._interp = interpolate.interp1d(
            self._log_x,
            self._log_s,
            kind="cubic",
            fill_value="extrapolate",
        )

    def __call__(self, x: float | npt.NDArray[np.float64]) -> float | npt.NDArray[np.float64]:
        """Evaluate s(x).

        Parameters
        ----------
        x : float or ndarray
            Position(s) at which to evaluate s.

        Returns
        -------
        float or ndarray
            Value(s) of s(x).
        """
        x_arr = np.atleast_1d(np.asarray(x, dtype=np.float64))
        result = np.zeros_like(x_arr)

        for i, xi in enumerate(x_arr):
            if self.x_min <= xi <= self.x_max:
                # Use interpolation
                result[i] = 10 ** self._interp(np.log10(xi))
            else:
                # Compute directly
                result[i] = compute_s_direct(xi, self.gamma)

        if np.isscalar(x):
            return float(result[0])
        return result


class MInterpolator:
    """Interpolator for m(x) = 2*pi integral_0^x dx' s(x') x'.

    Uses the SInterpolator and exploits logarithmic grid for integration.

    Parameters
    ----------
    s_interp : SInterpolator
        Interpolator for s(x).
    """

    def __init__(self, s_interp: SInterpolator) -> None:
        """Initialize the interpolator for m(x)."""
        self.s_interp = s_interp

        # Use the same grid as s
        self._x_grid = s_interp._x_grid.copy()
        self._log_x = s_interp._log_x.copy()

        # Compute m on the grid using cumulative trapezoid in log space
        # m(x) = integral_0^x s(x') x' dx'
        # With log grid: dx = x * ln(10) * d(log x)
        # So: m(x) = ln(10) integral s(x') x'^2 d(log x')
        # Note: s(x) is already normalized so integral_0^inf s(x) x dx = 1

        s_vals = s_interp._s_grid
        integrand = s_vals * self._x_grid**2

        # Cumulative trapezoid integration in log space
        d_log_x = self._log_x[1] - self._log_x[0]
        cumulative = integrate.cumulative_trapezoid(integrand, dx=d_log_x, initial=0)
        self._m_grid = np.log(10) * cumulative

        # Build inverse interpolator (m -> x) for finding effective radius
        # Only use unique, monotonically increasing m values
        valid = np.diff(self._m_grid) > 0
        valid = np.concatenate([[True], valid])
        m_valid = self._m_grid[valid]
        x_valid = self._x_grid[valid]

        self._inverse_interp = interpolate.interp1d(
            m_valid,
            np.log10(x_valid),
            kind="cubic",
            fill_value="extrapolate",
        )

        # Forward interpolator
        self._forward_interp = interpolate.interp1d(
            self._log_x,
            self._m_grid,
            kind="cubic",
            fill_value="extrapolate",
        )

    def __call__(self, x: float | npt.NDArray[np.float64]) -> float | npt.NDArray[np.float64]:
        """Evaluate m(x).

        Parameters
        ----------
        x : float or ndarray
            Position(s) at which to evaluate m.

        Returns
        -------
        float or ndarray
            Value(s) of m(x).
        """
        x_arr = np.atleast_1d(x)
        result: npt.NDArray[np.float64] = np.asarray(
            self._forward_interp(np.log10(x_arr)),
            dtype=np.float64,
        )

        if np.isscalar(x):
            return float(result[0])
        return result

    def inverse(self, m: float) -> float:
        """Find x such that m(x) = m.

        Parameters
        ----------
        m : float
            Target mass fraction.

        Returns
        -------
        float
            Position x where m(x) = m.
        """
        return float(10 ** self._inverse_interp(m))


def compute_effective_radius(gamma: float, normalize: bool = True) -> float:
    """Compute the effective radius x_e where m(x_e) = 0.5.

    Parameters
    ----------
    gamma : float
        Profile parameter.

    Returns
    -------
    float
        Effective radius x_e.
    """
    s_interp = SInterpolator(gamma)
    m_interp = MInterpolator(s_interp)
    effective_radius = m_interp.inverse(0.5)

    if normalize:
        half_mass_radius = 1 / (2 ** (1 / (3 - gamma)) - 1)
        return effective_radius / half_mass_radius
    return effective_radius


@register_asset("effective_radius_gamma.png")
def plot_effective_radius_gamma(*, path: Path | None = None) -> None:
    """Plot effective radius as a function of gamma.

    Parameters
    ----------
    path : Path or None
        If provided, save the figure to this path. Otherwise display it.
    """
    # Use 4 points for now while testing
    gamma_min, gamma_max = 0.5, 2.7
    gamma_values = np.linspace(gamma_min, gamma_max, 30)
    xe_values = np.array([compute_effective_radius(g) for g in gamma_values])

    fig: Figure
    ax: Axes
    fig, ax = plt.subplots(figsize=(6.5, 4.5))

    ax.plot(gamma_values, xe_values, "-", color="#000000", linewidth=2, markersize=8)
    ax.hlines(0.75, xmin=0, xmax=gamma_max, colors="#666666", linestyles="--", linewidth=2.0)

    ax.set_xlabel(r"$\gamma$")
    ax.set_ylabel(r"$R_e/r_{1/2}$")
    ax.set_xlim(gamma_min, gamma_max)
    ax.set_ylim(0.5, 0.9)

    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.tick_params(which="minor", length=3, color="gray", direction="in")
    ax.tick_params(which="major", length=6, direction="in")
    ax.tick_params(top=True, right=True, which="both")

    ax.grid(True, alpha=0.3)

    plt.tight_layout()

    if path:
        save_figure_if_changed(fig, path, dpi=150, bbox_inches="tight")
        plt.close(fig)
    else:
        plt.show()
