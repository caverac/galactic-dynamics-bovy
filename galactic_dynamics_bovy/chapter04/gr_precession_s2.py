"""GR precession of star S2 via direct 1PN orbit integration.

Computes the GR apsidal precession of star S2 orbiting Sgr A* using
the full first post-Newtonian (1PN) equation of motion in Cartesian 2D,
then compares with the effective-potential quadrature approach from
Problem 4.9 and the analytical approximation.

Parameters:
    M_bullet = 4.25e6 Msun, e = 0.884649, a = 1031.31 AU
"""

from dataclasses import dataclass
from pathlib import Path
from typing import Protocol

from matplotlib import pyplot as plt
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.ticker import AutoMinorLocator
import numpy as np
import numpy.typing as npt
from rich.console import Console
from scipy.integrate import OdeSolution, solve_ivp

from galactic_dynamics_bovy.chapter04.gr_precession import compute_delta_psi
from galactic_dynamics_bovy.utils.assets import register_asset, save_figure_if_changed
from galactic_dynamics_bovy.utils.units import SolarUnits


class OdeResult(Protocol):
    """Structural type for the object returned by ``solve_ivp``."""

    t: npt.NDArray[np.float64]
    y: npt.NDArray[np.float64]
    t_events: list[npt.NDArray[np.float64]]
    sol: OdeSolution
    success: bool


# ---------------------------------------------------------------------------
# S2 physical constants (all SI: m, s)
# ---------------------------------------------------------------------------
M_BULLET = 4.25e6  # black hole mass in solar masses
E_S2 = 0.884649  # orbital eccentricity
A_S2 = 1031.31 * SolarUnits.AU  # semi-major axis in metres
GM_BH = M_BULLET * SolarUnits.GM  # gravitational parameter (m^3/s^2)
C = SolarUnits.c  # speed of light (m/s)


def _eom_1pn(
    _t: float,
    y: npt.NDArray[np.float64],
    GM: float,
    c: float,
) -> list[float]:
    """Right-hand side of the 1PN equation of motion in Cartesian 2D.

    State vector: y = [x, y, vx, vy].

    The 1PN acceleration is:
        a_i = -GM*x_i/r^3 + GM/(c^2*r^2) * [(4GM/r^2 - v^2/r)*x_i + 4*v_r*v_i]
    """
    x, y_coord, vx, vy = y
    r2 = x * x + y_coord * y_coord
    r = np.sqrt(r2)
    r3 = r * r2
    v2 = vx * vx + vy * vy
    v_r = (x * vx + y_coord * vy) / r

    c2 = c * c

    # Newtonian acceleration
    ax = -GM * x / r3
    ay = -GM * y_coord / r3

    # 1PN correction
    coeff = GM / (c2 * r2)
    bracket_r = 4.0 * GM / r2 - v2 / r
    ax += coeff * (bracket_r * x + 4.0 * v_r * vx)
    ay += coeff * (bracket_r * y_coord + 4.0 * v_r * vy)

    return [vx, vy, ax, ay]


@dataclass(frozen=True)
class _PericenterEvent:
    """Callable event for ``solve_ivp`` detecting pericenter passages.

    ``direction=+1`` triggers when v_r crosses zero from negative to positive,
    i.e. the star passes through pericenter.
    """

    terminal: bool = False
    direction: int = 1

    def __call__(
        self,
        _t: float,
        y: npt.NDArray[np.float64],
    ) -> float:
        x, y_coord, vx, vy = y
        return float(x * vx + y_coord * vy)


_pericenter_event = _PericenterEvent()


def integrate_orbit_1pn(
    *,
    GM: float = GM_BH,
    c: float = C,
    e: float = E_S2,
    a: float = A_S2,
    n_orbits: float = 1.5,
    rtol: float = 1e-12,
    atol: float = 1e-12,
) -> OdeResult:
    """Integrate the 1PN orbit of S2 using DOP853.

    Parameters
    ----------
    GM : float
        Gravitational parameter (m^3/s^2).
    c : float
        Speed of light (m/s).
    e : float
        Orbital eccentricity.
    a : float
        Semi-major axis (m).
    n_orbits : float
        Number of Kepler periods to integrate.
    rtol, atol : float
        Relative and absolute tolerances for the integrator.

    Returns
    -------
    OdeResult
        Solution object from solve_ivp, with t_events populated.
    """
    # Initial conditions at pericenter
    r_peri = a * (1.0 - e)
    v_peri = np.sqrt(GM / a * (1.0 + e) / (1.0 - e))

    y0 = [r_peri, 0.0, 0.0, v_peri]

    # Kepler period for time span
    T_kepler = 2.0 * np.pi * np.sqrt(a**3 / GM)
    t_span = (0.0, n_orbits * T_kepler)

    # Closure captures GM and c, giving solve_ivp the (t, y) -> ... signature
    # it expects (avoids the poorly-typed ``args`` parameter).
    def eom(t: float, y: npt.NDArray[np.float64]) -> list[float]:
        return _eom_1pn(t, y, GM, c)

    sol = solve_ivp(
        eom,
        t_span,
        y0,
        method="DOP853",
        rtol=rtol,
        atol=atol,
        events=_pericenter_event,
        dense_output=True,
        max_step=T_kepler / 500.0,
    )
    return sol  # type: ignore[return-value]


def compute_precession_1pn(
    *,
    GM: float = GM_BH,
    c: float = C,
    e: float = E_S2,
    a: float = A_S2,
) -> float:
    """Measure GR precession per orbit from a direct 1PN integration.

    Integrates the orbit for ~1.5 periods, detects the next pericenter via
    the event function, and measures the angle of the position vector there.

    Returns
    -------
    float
        Precession per orbit (rad), i.e. Delta-psi - 2*pi.
    """
    sol = integrate_orbit_1pn(GM=GM, c=c, e=e, a=a, n_orbits=1.5)

    # First event after t=0 is the next pericenter
    t_events = sol.t_events[0]
    # Skip events very close to t=0 (the initial pericenter)
    T_kepler = 2.0 * np.pi * np.sqrt(a**3 / GM)
    mask = t_events > 0.1 * T_kepler
    t_peri = t_events[mask][0]

    state = sol.sol(t_peri)
    x_peri, y_peri = state[0], state[1]
    angle = np.arctan2(y_peri, x_peri)

    # The precession is the angle at the next pericenter
    # (which should be slightly > 0 for prograde precession)
    if angle < 0:
        angle += 2.0 * np.pi
    return float(angle)


def _setup_axes(ax: Axes) -> None:
    """Apply standard formatting to an axes object."""
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.tick_params(which="minor", length=3, color="gray", direction="in")
    ax.tick_params(which="major", length=6, direction="in")
    ax.tick_params(top=True, right=True, which="both")
    ax.grid(True, alpha=0.3)


@register_asset("p04_11_gr_precession_s2.png")
def plot_gr_precession_s2(path: Path | None = None) -> None:
    """Plot r(t) for the 1PN and Newtonian S2 orbits over several periods.

    The growing phase offset between the two curves reveals the GR
    apsidal precession.  Pericenter passages of the 1PN orbit are marked.

    Parameters
    ----------
    path : Path or None
        If provided, save the figure to this path. Otherwise display it.
    """
    n_orb = 5.5
    sol_1pn = integrate_orbit_1pn(n_orbits=n_orb)
    sol_newton = integrate_orbit_1pn(n_orbits=n_orb, c=1e30)

    T_kepler = 2.0 * np.pi * np.sqrt(A_S2**3 / GM_BH)
    sec_per_yr = 365.25 * 86400.0

    # 1PN trajectory
    t_yr = sol_1pn.t / sec_per_yr
    x = sol_1pn.y[0]
    y = sol_1pn.y[1]
    r_au = np.sqrt(x**2 + y**2) / SolarUnits.AU

    # Newtonian trajectory
    t_yr_n = sol_newton.t / sec_per_yr
    x_n = sol_newton.y[0]
    y_n = sol_newton.y[1]
    r_au_n = np.sqrt(x_n**2 + y_n**2) / SolarUnits.AU

    fig: Figure
    ax: Axes
    fig, ax = plt.subplots(figsize=(6.0, 4.0))
    fig.subplots_adjust(left=0.10, right=0.98, top=0.96, bottom=0.14)

    ax.plot(t_yr_n, r_au_n, "-", color="gray", linewidth=1, label="Newtonian")
    ax.plot(t_yr, r_au, "-", color="black", linewidth=2, label="1PN")

    # Mark 1PN pericenter passages
    t_events = sol_1pn.t_events[0]
    mask = t_events > 0.1 * T_kepler
    for t_peri in t_events[mask]:
        state = sol_1pn.sol(t_peri)
        r_peri = np.sqrt(state[0] ** 2 + state[1] ** 2) / SolarUnits.AU
        ax.plot(
            t_peri / sec_per_yr,
            r_peri,
            "o",
            color="black",
            markersize=4,
            zorder=5,
        )

    ax.set_xlabel(r"$t$ (yr)")
    ax.set_ylabel(r"$r$ (AU)")
    ax.legend(frameon=False, loc="upper right")
    _setup_axes(ax)

    if path:
        save_figure_if_changed(fig, path, dpi=150, bbox_inches="tight")
        plt.close(fig)
    else:
        plt.show()


if __name__ == "__main__":  # pragma: no cover
    console = Console()
    console.print("[bold]Problem 4.11: GR precession of S2 via 1PN integration[/bold]\n")

    # 1. Direct 1PN integration
    prec_1pn = compute_precession_1pn()
    prec_1pn_arcmin = np.degrees(prec_1pn) * 60.0

    # 2. Effective potential quadrature (from Problem 4.9)
    r_p = A_S2 * (1.0 - E_S2)
    prec_quad = compute_delta_psi(r_p, E_S2, GM=GM_BH, c=C) - 2.0 * np.pi
    prec_quad_arcmin = np.degrees(prec_quad) * 60.0

    # 3. Analytical: 6*pi*GM / (c^2 * a * (1 - e^2))
    prec_analytical = 6.0 * np.pi * GM_BH / (C**2 * A_S2 * (1.0 - E_S2**2))
    prec_analytical_arcmin = np.degrees(prec_analytical) * 60.0

    console.print("[green]Three-way comparison of GR precession per orbit:[/green]\n")
    console.print(f"  1PN integration:        {prec_1pn_arcmin:10.4f} arcmin")
    console.print(f"  Effective potential:     {prec_quad_arcmin:10.4f} arcmin")
    console.print(f"  Analytical (1st order):  {prec_analytical_arcmin:10.4f} arcmin")
    console.print()
    console.print(f"  1PN vs analytical:      {abs(prec_1pn - prec_analytical) / prec_analytical * 100:.4f}%")
    console.print(f"  Quadrature vs analytical: {abs(prec_quad - prec_analytical) / prec_analytical * 100:.4f}%")

    plot_gr_precession_s2()
