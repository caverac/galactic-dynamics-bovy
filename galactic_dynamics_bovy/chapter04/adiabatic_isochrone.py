"""Adiabatic invariance of actions in an isochrone potential.

Demonstrates that when the mass of an isochrone potential is doubled
**slowly** compared to the orbital period, the radial action J_r is
conserved (adiabatic invariant). When the mass grows **fast**, J_r changes.

The potential amplitude ramps from 1 to 2 via a smooth tanh profile:

    A(t) = 1.5 + 0.5 * tanh((t - t_mid) / tau)

where tau controls the ramp speed relative to the orbital period.
"""

from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path
import warnings

from galpy.actionAngle import actionAngleIsochrone
from galpy.orbit import Orbit
from galpy.potential import IsochronePotential, TimeDependentAmplitudeWrapperPotential, vcirc
from matplotlib import pyplot as plt
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.ticker import AutoMinorLocator
import numpy as np
import numpy.typing as npt

from galactic_dynamics_bovy.utils.assets import register_asset, save_figure_if_changed


@dataclass
class AdiabaticResult:
    """Results from an adiabatic mass-ramp experiment.

    Attributes
    ----------
    t : ndarray
        Time array.
    R : ndarray
        Galactocentric radius at each time step.
    Jr : ndarray
        Radial action computed at each snapshot time.
    Lz : ndarray
        Angular momentum at each snapshot time.
    amplitude : ndarray
        Potential amplitude A(t) at each snapshot time.
    T_orbit : float
        Estimated orbital period in the initial potential.
    tau : float
        Ramp timescale parameter.
    """

    t: npt.NDArray[np.float64]
    R: npt.NDArray[np.float64]
    Jr: npt.NDArray[np.float64]
    Lz: npt.NDArray[np.float64]
    amplitude: npt.NDArray[np.float64]
    T_orbit: float
    tau: float


def make_amplitude_func(tau: float, t_mid: float) -> Callable[[float], float]:
    """Create a smooth tanh ramp function from 1 to 2.

    Parameters
    ----------
    tau : float
        Ramp timescale (controls transition speed).
    t_mid : float
        Midpoint time of the ramp.

    Returns
    -------
    Callable[[float], float]
        Amplitude function A(t) = 1.5 + 0.5 * tanh((t - t_mid) / tau).
    """

    def amplitude(t: float) -> float:
        return float(1.5 + 0.5 * np.tanh((t - t_mid) / tau))

    return amplitude


def estimate_period(b: float = 1.0) -> float:
    """Estimate the orbital period for a mildly eccentric orbit.

    Integrates a short orbit in a static isochrone potential and returns
    the radial period via galpy's ``Orbit.Tp``.

    Parameters
    ----------
    b : float, optional
        Isochrone scale parameter. Default is 1.0.

    Returns
    -------
    float
        Estimated radial period in natural units.
    """
    base_pot = IsochronePotential(b=b, amp=1.0, normalize=False)
    vc = float(vcirc(base_pot, 1.0))
    o = Orbit([1.0, 0.05, 0.9 * vc, 0.0])
    t_short = np.linspace(0.0, 50.0, 5001)
    o.integrate(t_short, base_pot)
    with warnings.catch_warnings(), np.errstate(invalid="ignore"):
        warnings.simplefilter("ignore", RuntimeWarning)
        return float(o.Tp(pot=base_pot))


def run_adiabatic_experiment(
    tau: float,
    *,
    b: float = 1.0,
    n_orbits: int = 40,
    n_steps: int = 10001,
    n_snapshots: int = 200,
) -> AdiabaticResult:
    """Run an orbit in an isochrone potential whose mass ramps from 1 to 2.

    Parameters
    ----------
    tau : float
        Ramp timescale for the tanh amplitude function.
    b : float, optional
        Isochrone scale parameter. Default is 1.0.
    n_orbits : int, optional
        Number of orbital periods to integrate. Default is 40.
    n_steps : int, optional
        Number of integration time steps. Default is 10001.
    n_snapshots : int, optional
        Number of snapshot times for action evaluation. Default is 200.

    Returns
    -------
    AdiabaticResult
        Orbit trajectory and action time series.
    """
    T_orbit = estimate_period(b=b)
    t_max = n_orbits * T_orbit
    t_mid = t_max / 2.0

    base_pot = IsochronePotential(b=b, amp=1.0, normalize=False)
    amp_func = make_amplitude_func(tau, t_mid)
    td_pot = TimeDependentAmplitudeWrapperPotential(
        amp=1.0,
        A=amp_func,
        pot=base_pot,
    )

    vc = float(vcirc(base_pot, 1.0))
    o = Orbit([1.0, 0.05, 0.9 * vc, 0.0])
    t = np.linspace(0.0, t_max, n_steps)
    o.integrate(t, td_pot)

    R_full = o.R(t)

    # Compute actions at snapshot times
    t_snap = np.linspace(0.0, t_max, n_snapshots)
    Jr_arr = np.empty(n_snapshots)
    Lz_arr = np.empty(n_snapshots)
    amp_arr = np.empty(n_snapshots)

    # galpy's actionAngleSpherical emits RuntimeWarnings for planar orbits
    # (z=0) due to division by zero in inclination angle computation.
    # These are harmless and expected for our 2D orbit setup.
    with warnings.catch_warnings(), np.errstate(invalid="ignore"):
        warnings.simplefilter("ignore", RuntimeWarning)
        for i, ti in enumerate(t_snap):
            Ri = float(o.R(ti))
            vRi = float(o.vR(ti))
            vTi = float(o.vT(ti))
            amp_i = amp_func(ti)
            amp_arr[i] = amp_i

            snapshot_pot = IsochronePotential(b=b, amp=amp_i, normalize=False)
            aai = actionAngleIsochrone(ip=snapshot_pot)
            jr, lz, _ = aai(Ri, vRi, vTi, 0.0, 0.0)
            Jr_arr[i] = np.asarray(jr).flat[0]
            Lz_arr[i] = np.asarray(lz).flat[0]

    return AdiabaticResult(
        t=t,
        R=R_full,
        Jr=Jr_arr,
        Lz=Lz_arr,
        amplitude=amp_arr,
        T_orbit=T_orbit,
        tau=tau,
    )


def _setup_axes(ax: Axes) -> None:
    """Apply standard formatting to an axes object."""
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.tick_params(which="minor", length=3, color="gray", direction="in")
    ax.tick_params(which="major", length=6, direction="in")
    ax.tick_params(top=True, right=True, which="both")
    ax.grid(True, alpha=0.3)


@register_asset("p04_07_adiabatic_isochrone.png")
def plot_adiabatic_comparison(path: Path | None = None) -> None:
    """Plot 2x2 comparison of fast vs slow adiabatic mass ramp.

    Top row: fast ramp (tau ~ 0.1 * T_orbit).
    Bottom row: slow ramp (tau ~ 10 * T_orbit).
    Left column: R(t). Right column: J_r(t).

    Parameters
    ----------
    path : Path or None
        If provided, save the figure to this path. Otherwise display it.
    """
    T_orbit = estimate_period()
    tau_fast = 0.1 * T_orbit
    tau_slow = 10.0 * T_orbit

    result_fast = run_adiabatic_experiment(tau_fast)
    result_slow = run_adiabatic_experiment(tau_slow)

    fig: Figure
    axes: npt.NDArray[np.object_]
    fig, axes = plt.subplots(
        2,
        2,
        figsize=(7.0, 5.0),
        sharex=True,
    )
    fig.subplots_adjust(hspace=0.08, wspace=0.08, left=0.10, right=0.98, top=0.98, bottom=0.08)
    ax_jr_slow: Axes = axes[0, 0]
    ax_jr_fast: Axes = axes[0, 1]
    ax_r_slow: Axes = axes[1, 0]
    ax_r_fast: Axes = axes[1, 1]

    tau_label = r"$\tau = {:.1f}\,T_{{\rm orb}}$"

    # Bottom row: R(t)
    for ax, result, is_left in [
        (ax_r_slow, result_slow, True),
        (ax_r_fast, result_fast, False),
    ]:
        t_T = result.t / T_orbit
        t_mid = t_T[-1] / 2.0
        ramp_lo = t_mid - 2.0 * result.tau / T_orbit
        ramp_hi = t_mid + 2.0 * result.tau / T_orbit
        ax.plot(t_T, result.R, "-", color="#000000", linewidth=0.5)
        ax.axvspan(ramp_lo, ramp_hi, alpha=0.15, color="#4477AA")
        ax.set_xlabel(r"$t\,/\,T_{\rm orb}$")
        if is_left:
            ax.set_ylabel(r"$R\,/\,b$")
        else:
            ax.tick_params(labelleft=False)
        ax.text(
            0.97,
            0.95,
            tau_label.format(result.tau / T_orbit),
            transform=ax.transAxes,
            ha="right",
            va="top",
            fontsize=9,
        )
        _setup_axes(ax)

    # Top row: Jr(t)
    for ax, result, is_left in [
        (ax_jr_slow, result_slow, True),
        (ax_jr_fast, result_fast, False),
    ]:
        n_snap = len(result.Jr)
        t_snap = np.linspace(result.t[0], result.t[-1], n_snap) / T_orbit
        t_mid = t_snap[-1] / 2.0
        ramp_lo = t_mid - 2.0 * result.tau / T_orbit
        ramp_hi = t_mid + 2.0 * result.tau / T_orbit
        ax.plot(t_snap, result.Jr, "-", color="#000000", linewidth=1.0)
        ax.axhline(result.Jr[0], linestyle="--", color="#737373", linewidth=0.8, label=r"$J_r(0)$")
        ax.axvspan(ramp_lo, ramp_hi, alpha=0.15, color="#4477AA")
        if is_left:
            ax.set_ylabel(r"$J_r\,/\,(GM_0b)^{1/2}$")
        else:
            ax.tick_params(labelleft=False)
        ax.tick_params(labelbottom=False)
        ax.text(
            0.97,
            0.95,
            tau_label.format(result.tau / T_orbit),
            transform=ax.transAxes,
            ha="right",
            va="top",
            fontsize=9,
        )
        _setup_axes(ax)

    # Unify y-ranges within each row
    r_lo = min(ax_r_slow.get_ylim()[0], ax_r_fast.get_ylim()[0])
    r_hi = max(ax_r_slow.get_ylim()[1], ax_r_fast.get_ylim()[1])
    ax_r_slow.set_ylim(r_lo, r_hi)
    ax_r_fast.set_ylim(r_lo, r_hi)

    jr_lo = min(ax_jr_slow.get_ylim()[0], ax_jr_fast.get_ylim()[0])
    jr_hi = max(ax_jr_slow.get_ylim()[1], ax_jr_fast.get_ylim()[1])
    ax_jr_slow.set_ylim(jr_lo, jr_hi)
    ax_jr_fast.set_ylim(jr_lo, jr_hi)

    if path:
        save_figure_if_changed(fig, path, dpi=150, bbox_inches="tight")
        plt.close(fig)
    else:
        plt.show()


if __name__ == "__main__":  # pragma: no cover
    plot_adiabatic_comparison()
