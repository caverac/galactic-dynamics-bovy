"""NFW concentration from enclosed mass and escape velocity constraints.

This module solves for the NFW concentration parameter c given constraints
on enclosed mass and/or escape velocity (Problem 3.2).

Two approaches are provided:

1. **Using Dutton & Maccio (2014) relation** (Equation 2.66):
   Given M(<r), solve for c using the concentration-mass relation:
       log c = 0.905 - 0.101 log(M_200 / 10^12 h^-1 M_sun)

2. **Using escape velocity as constraint**:
   Given both M(<r) and v_esc(r), solve for (c, M_200) independently.
   This uses two observational constraints instead of a theoretical relation.

The enclosed mass for NFW is:
    M(<r) = M_200 * f(r/a) / f(c)

The NFW potential is:
    Phi(r) = -G M_200 / (r f(c)) * ln(1 + r/a)

The escape velocity to R_200 is:
    v_esc^2(r) = 2[Phi(R_200) - Phi(r)]

where f(x) = ln(1+x) - x/(1+x) and a = R_200/c is the scale radius.
"""

from dataclasses import dataclass

import numpy as np
import numpy.typing as npt
from scipy import optimize

from galactic_dynamics_bovy.utils.units import GalacticUnits


def nfw_f(x: float) -> float:
    """Compute the NFW mass function f(x) = ln(1+x) - x/(1+x).

    Parameters
    ----------
    x : float
        Dimensionless radius (r/a or concentration c).

    Returns
    -------
    float
        Value of f(x).
    """
    return float(np.log(1 + x) - x / (1 + x))


def m200_from_concentration(c: float, h: float = 0.67) -> float:
    """Compute M_200 from concentration using Dutton & Maccio (2014).

    Inverts Equation 2.66:
        log c = 0.905 - 0.101 log(M_200 / 10^12 h^-1 M_sun)

    Parameters
    ----------
    c : float
        Concentration parameter.
    h : float, optional
        Dimensionless Hubble parameter. Default is 0.67.

    Returns
    -------
    float
        M_200 in units of 10^10 M_sun.
    """
    # log(M_200 / 10^12 h^-1 M_sun) = (0.905 - log(c)) / 0.101
    log_m_ratio = (0.905 - np.log10(c)) / 0.101
    # M_200 = 10^12 h^-1 M_sun * 10^log_m_ratio
    # In units of 10^10 M_sun: M_200 = 100 * h^-1 * 10^log_m_ratio
    return float(100.0 / h * 10**log_m_ratio)


def r200_from_m200(m200: float, rho_crit: float) -> float:
    """Compute R_200 from M_200.

    Uses the definition M_200 = (4pi/3) * 200 * rho_crit * R_200^3.

    Parameters
    ----------
    m200 : float
        Virial mass in units of 10^10 M_sun.
    rho_crit : float
        Critical density in units of 10^10 M_sun / kpc^3.

    Returns
    -------
    float
        R_200 in kpc.
    """
    return float((3 * m200 / (800 * np.pi * rho_crit)) ** (1.0 / 3.0))


def nfw_potential_interior(r: float, m200: float, c: float, r200: float, G: float) -> float:
    """Potential from mass inside r (Term 1).

    Phi_interior(r) = -G M(<r) / r

    The enclosed mass M(<r) acts as a point mass at the center.

    Parameters
    ----------
    r : float
        Radius in kpc.
    m200 : float
        Virial mass in units of 10^10 M_sun.
    c : float
        Concentration parameter.
    r200 : float
        Virial radius in kpc.
    G : float
        Gravitational constant in kpc (km/s)^2 / (10^10 M_sun).

    Returns
    -------
    float
        Interior potential contribution in (km/s)^2.
    """
    a = r200 / c
    fc = nfw_f(c)
    m_enclosed = m200 * nfw_f(r / a) / fc
    return float(-G * m_enclosed / r)


def nfw_potential_exterior(r: float, m200: float, c: float, r200: float, G: float) -> float:
    """Potential from mass between r and R_200 (Term 2).

    Phi_exterior(r) = -G M_200 / f(c) x [1/(a+r) - 1/(a+R_200)]

    This is the contribution from being inside the shell of mass
    between r and R_200. Derived from -4*pi*G integral_r^R_200 rho(r') r' dr'.

    Parameters
    ----------
    r : float
        Radius in kpc.
    m200 : float
        Virial mass in units of 10^10 M_sun.
    c : float
        Concentration parameter.
    r200 : float
        Virial radius in kpc.
    G : float
        Gravitational constant in kpc (km/s)^2 / (10^10 M_sun).

    Returns
    -------
    float
        Exterior potential contribution in (km/s)^2.
    """
    a = r200 / c
    fc = nfw_f(c)
    return float(-G * m200 / fc * (1 / (a + r) - 1 / (a + r200)))


def nfw_potential(r: float, m200: float, c: float, r200: float, G: float) -> float:
    """Compute NFW potential at radius r for halo truncated at R_200.

    Phi(r) = Phi_interior(r) + Phi_exterior(r)

    where:
        Phi_interior(r) = -G M(<r) / r           (mass inside r)
        Phi_exterior(r) = -G M_200/f(c) x [1/(a+r) - 1/(a+R_200)]  (shell r to R_200)

    Parameters
    ----------
    r : float
        Radius in kpc.
    m200 : float
        Virial mass in units of 10^10 M_sun.
    c : float
        Concentration parameter.
    r200 : float
        Virial radius in kpc.
    G : float
        Gravitational constant in kpc (km/s)^2 / (10^10 M_sun).

    Returns
    -------
    float
        Potential in (km/s)^2.
    """
    phi_int = nfw_potential_interior(r, m200, c, r200, G)
    phi_ext = nfw_potential_exterior(r, m200, c, r200, G)
    return float(phi_int + phi_ext)


def nfw_vesc_to_r200(r: float, m200: float, c: float, r200: float, G: float) -> float:
    """Compute escape velocity from r to R_200 for an NFW profile.

    v_esc^2(r) = 2[Phi(R_200) - Phi(r)]

    Parameters
    ----------
    r : float
        Radius in kpc.
    m200 : float
        Virial mass in units of 10^10 M_sun.
    c : float
        Concentration parameter.
    r200 : float
        Virial radius in kpc.
    G : float
        Gravitational constant in kpc (km/s)^2 / (10^10 M_sun).

    Returns
    -------
    float
        Escape velocity in km/s.
    """
    phi_r = nfw_potential(r, m200, c, r200, G)
    phi_r200 = nfw_potential(r200, m200, c, r200, G)
    vesc_squared = 2 * (phi_r200 - phi_r)
    return float(np.sqrt(max(vesc_squared, 0.0)))


def enclosed_mass_residual(
    c: float,
    r_inner: float,
    m_inner: float,
    rho_crit: float,
    h: float,
) -> float:
    """Compute residual for the enclosed mass constraint (Equation 3.2.2).

    The constraint is:
        M_inner = M_200(c) * f(r_inner * c / R_200) / f(c)

    Parameters
    ----------
    c : float
        Concentration parameter (to solve for).
    r_inner : float
        Inner radius in kpc.
    m_inner : float
        Enclosed mass within r_inner in units of 10^10 M_sun.
    rho_crit : float
        Critical density in units of 10^10 M_sun / kpc^3.
    h : float
        Dimensionless Hubble parameter.

    Returns
    -------
    float
        Residual (predicted M_inner - observed M_inner).
    """
    m200 = m200_from_concentration(c, h)
    r200 = r200_from_m200(m200, rho_crit)

    # x = r_inner / a = r_inner * c / R_200
    x = r_inner * c / r200

    # Predicted enclosed mass
    m_predicted = m200 * nfw_f(x) / nfw_f(c)

    return float(m_predicted - m_inner)


@dataclass
class SolvedNFWModel:
    """Dataclass to hold solved NFW concentration results."""

    c: float
    m200: float
    r200: float
    a: float
    v_esc_predicted: float


def solve_concentration(
    *,
    r_inner: float = 8.0,
    m_inner: float = 9.0,
    h: float = 0.67,
    c_bracket: tuple[float, float] = (5.0, 20.0),
) -> SolvedNFWModel:
    """Solve for NFW concentration given enclosed mass constraint.

    Solves Equation 3.2.2:
        M_inner = M_200(c) * f(r_inner * c / R_200) / f(c)

    Parameters
    ----------
    r_inner : float, optional
        Inner radius in kpc. Default is 8.0 (solar radius).
    m_inner : float, optional
        Enclosed mass within r_inner in units of 10^10 M_sun.
        Default is 9.0 (9 x 10^10 M_sun).
    h : float, optional
        Dimensionless Hubble parameter. Default is 0.67.
    c_bracket : tuple[float, float], optional
        Bracket for root finding. Default is (5.0, 20.0).

    Returns
    -------
    SolvedNFWModel
        - 'c': concentration parameter
        - 'm200': M_200 in 10^10 M_sun
        - 'r200': R_200 in kpc
        - 'a': scale radius in kpc
    """
    units = GalacticUnits(h=h)

    result = optimize.brentq(
        enclosed_mass_residual,
        c_bracket[0],
        c_bracket[1],
        args=(r_inner, m_inner, units.rho_crit, h),
    )
    c = float(result)

    m200 = m200_from_concentration(c, h)
    r200 = r200_from_m200(m200, units.rho_crit)
    a = r200 / c
    v_esc = nfw_vesc_to_r200(r_inner, m200, c, r200, units.G_kms)

    return SolvedNFWModel(c=c, m200=m200, r200=r200, a=a, v_esc_predicted=v_esc)


def _vesc_residuals(
    params: npt.NDArray[np.float64],
    r_inner: float,
    m_inner: float,
    v_esc: float,
    *,
    rho_crit: float,
    G: float,
) -> npt.NDArray[np.float64]:
    """Compute residuals for the two-constraint system.

    Parameters
    ----------
    params : ndarray
        Array [c, m200] - concentration and virial mass.
    r_inner : float
        Inner radius in kpc.
    m_inner : float
        Enclosed mass within r_inner in units of 10^10 M_sun.
    v_esc : float
        Escape velocity at r_inner in km/s.
    rho_crit : float
        Critical density in units of 10^10 M_sun / kpc^3.
    G : float
        Gravitational constant in kpc (km/s)^2 / (10^10 M_sun).

    Returns
    -------
    ndarray
        Array [mass_residual, vesc_residual].
    """
    c, m200 = float(params[0]), float(params[1])
    r200 = r200_from_m200(m200, rho_crit)
    a = r200 / c

    # Mass constraint: M(<r_inner) = m_inner
    x = r_inner / a
    m_predicted = m200 * nfw_f(x) / nfw_f(c)
    mass_residual = m_predicted - m_inner

    # Escape velocity constraint: v_esc(r_inner) = v_esc
    vesc_predicted = nfw_vesc_to_r200(r_inner, m200, c, r200, G)
    vesc_residual = vesc_predicted - v_esc

    return np.array([mass_residual, vesc_residual])


def solve_concentration_from_vesc(
    *,
    r_inner: float = 8.0,
    m_inner: float = 9.0,
    v_esc: float = 550.0,
    h: float = 0.67,
    initial_guess: tuple[float, float] = (10.0, 100.0),
) -> SolvedNFWModel:
    """Solve for NFW (c, M_200) using enclosed mass and escape velocity.

    Uses two observational constraints:
    1. M(<r_inner) = m_inner
    2. v_esc(r_inner) = v_esc (escape to R_200)

    This does NOT use the Dutton & Maccio concentration-mass relation.

    Parameters
    ----------
    r_inner : float, optional
        Inner radius in kpc. Default is 8.0 (solar radius).
    m_inner : float, optional
        Enclosed mass within r_inner in units of 10^10 M_sun.
        Default is 9.0 (9 x 10^10 M_sun).
    v_esc : float, optional
        Escape velocity at r_inner in km/s. Default is 550.0.
    h : float, optional
        Dimensionless Hubble parameter. Default is 0.67.
    initial_guess : tuple[float, float], optional
        Initial guess for (c, m200). Default is (10.0, 100.0).

    Returns
    -------
    SolvedNFWModel
        - 'c': concentration parameter
        - 'm200': M_200 in 10^10 M_sun
        - 'r200': R_200 in kpc
        - 'a': scale radius in kpc
        - 'v_esc_predicted': predicted escape velocity (should match v_esc)
    """
    units = GalacticUnits(h=h)

    def residuals(params: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
        return _vesc_residuals(params, r_inner, m_inner, v_esc, rho_crit=units.rho_crit, G=units.G_kms)

    result = optimize.fsolve(
        residuals,
        np.array(initial_guess),
        full_output=True,
    )
    solution, _, ier, mesg = result

    if ier != 1:
        raise RuntimeError(f"Failed to converge: {mesg}")

    c, m200 = float(solution[0]), float(solution[1])
    r200 = r200_from_m200(m200, units.rho_crit)
    a = r200 / c

    v_esc_predicted = nfw_vesc_to_r200(r_inner, m200, c, r200, units.G_kms)

    return SolvedNFWModel(c=c, m200=m200, r200=r200, a=a, v_esc_predicted=v_esc_predicted)
