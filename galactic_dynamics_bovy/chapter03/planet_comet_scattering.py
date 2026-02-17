"""Orbit scattering in an isochrone potential (Problem 4.3).

This module studies what happens to an orbit in an isochrone potential when
an instantaneous velocity kick is applied at pericenter. It investigates:

1. Whether the new pericenter is larger or smaller after the kick
2. The effect of radial-only vs tangential-only kicks
3. How the pericenter change depends on kick magnitude

The isochrone potential is defined as:
    Phi(r) = -G M / (b + sqrt(b^2 + r^2))

where b is the scale parameter.
"""

from dataclasses import dataclass
from pathlib import Path

from galpy.orbit import Orbit
from galpy.potential import IsochronePotential
from matplotlib import pyplot as plt
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.ticker import AutoMinorLocator
import numpy as np
import numpy.typing as npt
from scipy.signal import argrelmin

from galactic_dynamics_bovy.utils.assets import register_asset, save_figure_if_changed


@dataclass
class OrbitResult:
    """Bundle a galpy Orbit with extracted orbital quantities.

    Attributes
    ----------
    orbit : Orbit
        The integrated galpy Orbit object.
    r_peri : float
        Pericenter radius.
    r_apo : float
        Apocenter radius.
    energy : float
        Orbital energy.
    lz : float
        z-component of angular momentum.
    t : ndarray
        Time array used for integration.
    """

    orbit: Orbit
    r_peri: float
    r_apo: float
    energy: float
    lz: float
    t: npt.NDArray[np.float64]


@dataclass
class ScatteringResult:
    """Pair original and perturbed OrbitResult with kick parameters.

    Attributes
    ----------
    original : OrbitResult
        The unperturbed orbit.
    perturbed : OrbitResult
        The orbit after the velocity kick.
    delta_v_r : float
        Radial velocity kick applied.
    delta_v_t : float
        Tangential velocity kick applied.
    delta_r_peri : float
        Change in pericenter radius (perturbed - original).
    delta_r_peri_frac : float
        Fractional change in pericenter radius.
    """

    original: OrbitResult
    perturbed: OrbitResult
    delta_v_r: float
    delta_v_t: float
    delta_r_peri: float
    delta_r_peri_frac: float


def setup_isochrone_potential(*, b: float = 1.0) -> IsochronePotential:  # type: ignore[type-arg]
    """Create a galpy IsochronePotential with the given scale parameter.

    Parameters
    ----------
    b : float, optional
        Scale parameter of the isochrone potential. Default is 1.0.

    Returns
    -------
    IsochronePotential
        The configured potential in natural units.
    """
    return IsochronePotential(b=b, normalize=True)


def integrate_orbit(
    *,
    R: float,
    vR: float,
    vT: float,
    b: float = 1.0,
    t_max: float = 50.0,
    n_steps: int = 10001,
) -> OrbitResult:
    """Integrate an orbit in an isochrone potential.

    Parameters
    ----------
    R : float
        Initial galactocentric radius (natural units).
    vR : float
        Initial radial velocity (natural units).
    vT : float
        Initial tangential velocity (natural units).
    b : float, optional
        Isochrone scale parameter. Default is 1.0.
    t_max : float, optional
        Maximum integration time. Default is 50.0.
    n_steps : int, optional
        Number of time steps. Default is 10001.

    Returns
    -------
    OrbitResult
        Integrated orbit with extracted quantities.
    """
    pot = setup_isochrone_potential(b=b)
    t = np.linspace(0.0, t_max, n_steps)

    o = Orbit([R, vR, vT, 0.0])
    o.integrate(t, pot)

    r_values = o.R(t)
    r_peri = float(np.min(r_values))
    r_apo = float(np.max(r_values))
    energy = float(o.E(pot=pot))
    lz = float(o.Lz())

    return OrbitResult(orbit=o, r_peri=r_peri, r_apo=r_apo, energy=energy, lz=lz, t=t)


def find_pericenter_time(orbit_result: OrbitResult, *, occurrence: int = 1) -> float:
    """Find the time of the Nth pericenter passage.

    Detects local minima in r(t) using scipy.signal.argrelmin.

    Parameters
    ----------
    orbit_result : OrbitResult
        An integrated orbit result.
    occurrence : int, optional
        Which pericenter passage to find (1-indexed). Default is 1.

    Returns
    -------
    float
        Time of the requested pericenter passage.

    Raises
    ------
    ValueError
        If the requested occurrence is not found.
    """
    if occurrence < 1:
        raise ValueError(f"occurrence must be >= 1, got {occurrence}")

    r_values = orbit_result.orbit.R(orbit_result.t)
    minima_indices = argrelmin(r_values, order=5)[0]

    if len(minima_indices) < occurrence:
        raise ValueError(f"Only {len(minima_indices)} pericenter(s) found, " f"requested occurrence={occurrence}")

    idx = minima_indices[occurrence - 1]
    return float(orbit_result.t[idx])


def apply_kick_at_pericenter(
    orbit_result: OrbitResult,
    *,
    delta_v_r: float = 0.0,
    delta_v_t: float = 0.0,
    b: float = 1.0,
    t_max: float = 50.0,
    n_steps: int = 10001,
    pericenter_occurrence: int = 1,
) -> ScatteringResult:
    """Apply a velocity kick at pericenter and integrate the new orbit.

    Parameters
    ----------
    orbit_result : OrbitResult
        The original integrated orbit.
    delta_v_r : float, optional
        Radial velocity kick. Default is 0.0.
    delta_v_t : float, optional
        Tangential velocity kick. Default is 0.0.
    b : float, optional
        Isochrone scale parameter. Default is 1.0.
    t_max : float, optional
        Maximum integration time for the new orbit. Default is 50.0.
    n_steps : int, optional
        Number of time steps. Default is 10001.
    pericenter_occurrence : int, optional
        Which pericenter passage to apply the kick at. Default is 1.

    Returns
    -------
    ScatteringResult
        Comparison of original and perturbed orbits.
    """
    t_peri = find_pericenter_time(orbit_result, occurrence=pericenter_occurrence)

    o = orbit_result.orbit
    R_peri = float(o.R(t_peri))
    vR_peri = float(o.vR(t_peri))
    vT_peri = float(o.vT(t_peri))

    new_vR = vR_peri + delta_v_r
    new_vT = vT_peri + delta_v_t

    perturbed = integrate_orbit(
        R=R_peri,
        vR=new_vR,
        vT=new_vT,
        b=b,
        t_max=t_max,
        n_steps=n_steps,
    )

    delta_r_peri = perturbed.r_peri - orbit_result.r_peri
    delta_r_peri_frac = delta_r_peri / orbit_result.r_peri

    return ScatteringResult(
        original=orbit_result,
        perturbed=perturbed,
        delta_v_r=delta_v_r,
        delta_v_t=delta_v_t,
        delta_r_peri=delta_r_peri,
        delta_r_peri_frac=delta_r_peri_frac,
    )


@dataclass
class ParameterStudyResult:
    """Results from a parameter study sweep.

    Attributes
    ----------
    kicks : ndarray
        Array of kick magnitudes used.
    delta_peri_radial : list of ndarray
        Log pericenter ratio for radial kicks, one array per orbit.
    delta_peri_tangential : list of ndarray
        Log pericenter ratio for tangential kicks, one array per orbit.
    delta_apo_radial : list of ndarray
        Log apocenter ratio for radial kicks, one array per orbit.
    delta_apo_tangential : list of ndarray
        Log apocenter ratio for tangential kicks, one array per orbit.
    delta_ecc_radial : list of ndarray
        Log eccentricity ratio for radial kicks, one array per orbit.
    delta_ecc_tangential : list of ndarray
        Log eccentricity ratio for tangential kicks, one array per orbit.
    """

    kicks: npt.NDArray[np.float64]
    delta_peri_radial: list[npt.NDArray[np.float64]]
    delta_peri_tangential: list[npt.NDArray[np.float64]]
    delta_apo_radial: list[npt.NDArray[np.float64]]
    delta_apo_tangential: list[npt.NDArray[np.float64]]
    delta_ecc_radial: list[npt.NDArray[np.float64]]
    delta_ecc_tangential: list[npt.NDArray[np.float64]]


def _eccentricity(r_peri: float, r_apo: float) -> float:
    """Compute orbital eccentricity from pericenter and apocenter."""
    return (r_apo - r_peri) / (r_apo + r_peri)


def parameter_study(
    *,
    kick_magnitudes: npt.NDArray[np.float64] | None = None,
    initial_conditions: list[tuple[float, float, float]] | None = None,
    b: float = 1.0,
) -> ParameterStudyResult:
    """Sweep over kick magnitudes for multiple initial orbits.

    Parameters
    ----------
    kick_magnitudes : ndarray or None, optional
        Array of kick magnitudes to test. Default is np.linspace(0.01, 0.3, 30).
    initial_conditions : list of tuples or None, optional
        List of (R, vR, vT) initial conditions. Default provides two orbits.
    b : float, optional
        Isochrone scale parameter. Default is 1.0.

    Returns
    -------
    ParameterStudyResult
        Results containing log ratios for pericenter, apocenter, and eccentricity.
    """
    if kick_magnitudes is None:
        kick_magnitudes = np.linspace(0.01, 0.3, 30)

    if initial_conditions is None:
        initial_conditions = [
            (1.0, 0.3, 0.5),
            (0.5, 0.0, 0.8),
        ]

    delta_peri_radial: list[npt.NDArray[np.float64]] = []
    delta_peri_tangential: list[npt.NDArray[np.float64]] = []
    delta_apo_radial: list[npt.NDArray[np.float64]] = []
    delta_apo_tangential: list[npt.NDArray[np.float64]] = []
    delta_ecc_radial: list[npt.NDArray[np.float64]] = []
    delta_ecc_tangential: list[npt.NDArray[np.float64]] = []

    for R, vR, vT in initial_conditions:
        orbit = integrate_orbit(R=R, vR=vR, vT=vT, b=b)
        ecc_orig = _eccentricity(orbit.r_peri, orbit.r_apo)

        peri_r = np.empty(len(kick_magnitudes))
        peri_t = np.empty(len(kick_magnitudes))
        apo_r = np.empty(len(kick_magnitudes))
        apo_t = np.empty(len(kick_magnitudes))
        ecc_r = np.empty(len(kick_magnitudes))
        ecc_t = np.empty(len(kick_magnitudes))

        for i, dv in enumerate(kick_magnitudes):
            result_r = apply_kick_at_pericenter(orbit, delta_v_r=dv, b=b)
            peri_r[i] = np.log(result_r.perturbed.r_peri / orbit.r_peri)
            apo_r[i] = np.log(result_r.perturbed.r_apo / orbit.r_apo)
            ecc_r[i] = np.log(_eccentricity(result_r.perturbed.r_peri, result_r.perturbed.r_apo) / ecc_orig)

            result_t = apply_kick_at_pericenter(orbit, delta_v_t=dv, b=b)
            peri_t[i] = np.log(result_t.perturbed.r_peri / orbit.r_peri)
            apo_t[i] = np.log(result_t.perturbed.r_apo / orbit.r_apo)
            ecc_t[i] = np.log(_eccentricity(result_t.perturbed.r_peri, result_t.perturbed.r_apo) / ecc_orig)

        delta_peri_radial.append(peri_r)
        delta_peri_tangential.append(peri_t)
        delta_apo_radial.append(apo_r)
        delta_apo_tangential.append(apo_t)
        delta_ecc_radial.append(ecc_r)
        delta_ecc_tangential.append(ecc_t)

    return ParameterStudyResult(
        kicks=kick_magnitudes,
        delta_peri_radial=delta_peri_radial,
        delta_peri_tangential=delta_peri_tangential,
        delta_apo_radial=delta_apo_radial,
        delta_apo_tangential=delta_apo_tangential,
        delta_ecc_radial=delta_ecc_radial,
        delta_ecc_tangential=delta_ecc_tangential,
    )


def _setup_axes(ax: Axes) -> None:
    """Apply standard formatting to an axes object."""
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.tick_params(which="minor", length=3, color="gray", direction="in")
    ax.tick_params(which="major", length=6, direction="in")
    ax.tick_params(top=True, right=True, which="both")
    ax.grid(True, alpha=0.3)


@register_asset("p04_10_orbit_comparison.png")
def plot_orbit_comparison(path: Path | None = None) -> None:
    """Plot original vs kicked orbit: x-y traces and R(t).

    Creates a 1x2 figure comparing the original orbit with an orbit
    that received a radial velocity kick at pericenter.

    Parameters
    ----------
    path : Path or None
        If provided, save the figure to this path. Otherwise display it.
    """
    orbit = integrate_orbit(R=1.0, vR=0.3, vT=0.5, t_max=20.0)
    result = apply_kick_at_pericenter(orbit, delta_v_r=0.15, t_max=20.0)

    fig: Figure
    axes: npt.NDArray[np.object_]
    fig, axes = plt.subplots(1, 2, figsize=(6.5, 4.5))
    ax_xy: Axes = axes[0]
    ax_rt: Axes = axes[1]

    o_orig = result.original.orbit
    o_pert = result.perturbed.orbit

    # Show original orbit from pericenter onward so both start at same R
    t_peri = find_pericenter_time(result.original)
    t_full = result.original.t
    mask = t_full >= t_peri
    t_orig = t_full[mask] - t_peri  # shift so both start at t=0
    t_pert = result.perturbed.t

    # x-y orbit traces
    ax_xy.plot(
        o_orig.R(t_full[mask]) * np.cos(o_orig.phi(t_full[mask])),
        o_orig.R(t_full[mask]) * np.sin(o_orig.phi(t_full[mask])),
        "-",
        color="#000000",
        linewidth=1.5,
        label="Original",
    )
    ax_xy.plot(
        o_pert.R(t_pert) * np.cos(o_pert.phi(t_pert)),
        o_pert.R(t_pert) * np.sin(o_pert.phi(t_pert)),
        "--",
        color="#737373",
        linewidth=1.5,
        label="Radially kicked",
    )
    ax_xy.set_xlabel(r"$x$")
    ax_xy.set_ylabel(r"$y$")
    ax_xy.set_aspect("equal", adjustable="datalim")
    ax_xy.legend(frameon=False)
    _setup_axes(ax_xy)

    # R(t)
    ax_rt.plot(t_orig, o_orig.R(t_full[mask]), "-", color="#000000", linewidth=1.5)
    ax_rt.plot(t_pert, o_pert.R(t_pert), "--", color="#737373", linewidth=1.5)
    ax_rt.set_xlabel(r"$t$")
    ax_rt.set_ylabel(r"$R$")
    _setup_axes(ax_rt)

    plt.tight_layout()

    if path:
        save_figure_if_changed(fig, path, dpi=150, bbox_inches="tight")
        plt.close(fig)
    else:
        plt.show()


@register_asset("p04_10_parameter_study.png")
def plot_parameter_study(path: Path | None = None) -> None:
    """Plot parameter study comparing multiple initial orbits.

    1x2 figure: radial kicks (left) and tangential kicks (right),
    each showing results for different initial conditions.

    Parameters
    ----------
    path : Path or None
        If provided, save the figure to this path. Otherwise display it.
    """
    result = parameter_study()
    kicks = result.kicks
    delta_radial = result.delta_peri_radial
    delta_tangential = result.delta_peri_tangential

    orbit_labels = [r"$R=1.0,\, v_R=0.3,\, v_T=0.5$", r"$R=0.5,\, v_R=0.0,\, v_T=0.8$"]
    colors = ["#000000", "#737373"]
    styles = ["-", "--"]

    fig: Figure
    axes: npt.NDArray[np.object_]
    fig, axes = plt.subplots(1, 2, figsize=(6.5, 4.5), sharey=True)
    fig.subplots_adjust(wspace=0.08)
    ax_radial: Axes = axes[0]
    ax_tangential: Axes = axes[1]

    for i, label in enumerate(orbit_labels):
        ax_radial.plot(
            kicks,
            delta_radial[i],
            styles[i],
            color=colors[i],
            linewidth=1.5,
            label=label,
        )
        ax_tangential.plot(
            kicks,
            delta_tangential[i],
            styles[i],
            color=colors[i],
            linewidth=1.5,
            label=label,
        )

    ax_radial.set_xlabel(r"$|\Delta v_r|$")
    ax_radial.set_ylabel(r"$\mathrm{ln}\, r_p^{(k)} / r_p$")
    ax_radial.set_title("Radial kicks")
    ax_radial.axhline(0, color="gray", linewidth=0.5)
    ax_radial.legend(frameon=False)
    ax_radial.set_ylim(-0.12, 0.01)
    _setup_axes(ax_radial)

    ax_tangential.set_xlabel(r"$|\Delta v_T|$")
    ax_tangential.set_ylabel(r"$\mathrm{ln}\, r_p^{(k)} / r_p$")
    ax_tangential.set_title("Tangential kicks")
    ax_tangential.axhline(0, color="gray", linewidth=0.5)
    _setup_axes(ax_tangential)

    plt.tight_layout()

    if path:
        save_figure_if_changed(fig, path, dpi=150, bbox_inches="tight")
        plt.close(fig)
    else:
        plt.show()


@register_asset("p04_10_apocenter_study.png")
def plot_apocenter_study(path: Path | None = None) -> None:
    """Plot apocenter change vs kick magnitude.

    1x2 figure: radial kicks (left) and tangential kicks (right),
    each showing results for different initial conditions.

    Parameters
    ----------
    path : Path or None
        If provided, save the figure to this path. Otherwise display it.
    """
    result = parameter_study()
    kicks = result.kicks
    delta_radial = result.delta_apo_radial
    delta_tangential = result.delta_apo_tangential

    orbit_labels = [r"$R=1.0,\, v_R=0.3,\, v_T=0.5$", r"$R=0.5,\, v_R=0.0,\, v_T=0.8$"]
    colors = ["#000000", "#737373"]
    styles = ["-", "--"]

    fig: Figure
    axes: npt.NDArray[np.object_]
    fig, axes = plt.subplots(1, 2, figsize=(6.5, 4.5), sharey=True)
    fig.subplots_adjust(wspace=0.08)
    ax_radial: Axes = axes[0]
    ax_tangential: Axes = axes[1]

    for i, label in enumerate(orbit_labels):
        ax_radial.plot(
            kicks,
            delta_radial[i],
            styles[i],
            color=colors[i],
            linewidth=1.5,
            label=label,
        )
        ax_tangential.plot(
            kicks,
            delta_tangential[i],
            styles[i],
            color=colors[i],
            linewidth=1.5,
            label=label,
        )

    ax_radial.set_xlabel(r"$|\Delta v_r|$")
    ax_radial.set_ylabel(r"$\mathrm{ln}\, r_a^{(k)} / r_a$")
    ax_radial.set_title("Radial kicks")
    ax_radial.axhline(0, color="gray", linewidth=0.5)
    ax_radial.legend(frameon=False)
    _setup_axes(ax_radial)

    ax_tangential.set_xlabel(r"$|\Delta v_T|$")
    ax_tangential.set_ylabel(r"$\mathrm{ln}\, r_a^{(k)} / r_a$")
    ax_tangential.set_title("Tangential kicks")
    ax_tangential.axhline(0, color="gray", linewidth=0.5)
    _setup_axes(ax_tangential)

    plt.tight_layout()

    if path:
        save_figure_if_changed(fig, path, dpi=150, bbox_inches="tight")
        plt.close(fig)
    else:
        plt.show()


@register_asset("p04_10_eccentricity_study.png")
def plot_eccentricity_study(path: Path | None = None) -> None:
    """Plot eccentricity change vs kick magnitude.

    1x2 figure: radial kicks (left) and tangential kicks (right),
    each showing results for different initial conditions.

    Parameters
    ----------
    path : Path or None
        If provided, save the figure to this path. Otherwise display it.
    """
    result = parameter_study()
    kicks = result.kicks
    delta_radial = result.delta_ecc_radial
    delta_tangential = result.delta_ecc_tangential

    orbit_labels = [r"$R=1.0,\, v_R=0.3,\, v_T=0.5$", r"$R=0.5,\, v_R=0.0,\, v_T=0.8$"]
    colors = ["#000000", "#737373"]
    styles = ["-", "--"]

    fig: Figure
    axes: npt.NDArray[np.object_]
    fig, axes = plt.subplots(1, 2, figsize=(6.5, 4.5), sharey=True)
    fig.subplots_adjust(wspace=0.08)
    ax_radial: Axes = axes[0]
    ax_tangential: Axes = axes[1]

    for i, label in enumerate(orbit_labels):
        ax_radial.plot(
            kicks,
            delta_radial[i],
            styles[i],
            color=colors[i],
            linewidth=1.5,
            label=label,
        )
        ax_tangential.plot(
            kicks,
            delta_tangential[i],
            styles[i],
            color=colors[i],
            linewidth=1.5,
            label=label,
        )

    ax_radial.set_xlabel(r"$|\Delta v_r|$")
    ax_radial.set_ylabel(r"$\mathrm{ln}\, e^{(k)} / e$")
    ax_radial.set_title("Radial kicks")
    ax_radial.axhline(0, color="gray", linewidth=0.5)
    ax_radial.legend(frameon=False)
    _setup_axes(ax_radial)

    ax_tangential.set_xlabel(r"$|\Delta v_T|$")
    ax_tangential.set_ylabel(r"$\mathrm{ln}\, e^{(k)} / e$")
    ax_tangential.set_title("Tangential kicks")
    ax_tangential.axhline(0, color="gray", linewidth=0.5)
    _setup_axes(ax_tangential)

    plt.tight_layout()

    if path:
        save_figure_if_changed(fig, path, dpi=150, bbox_inches="tight")
        plt.close(fig)
    else:
        plt.show()


if __name__ == "__main__":  # pragma: no cover
    plot_orbit_comparison()
    plot_parameter_study()
    plot_apocenter_study()
    plot_eccentricity_study()
