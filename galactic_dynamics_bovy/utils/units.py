# pyright: reportAttributeAccessIssue=false
# pyright: reportOperatorIssue=false
# pylint: disable=no-member, invalid-name
# mypy: disable-error-code="attr-defined, union-attr, operator"
"""Physical constants and units for galactic dynamics.

This module provides a convenient interface to physical constants
in two unit systems commonly used in galactic dynamics:

- StellarUnits: (pc, Msun, Myr, km/s)
- GalacticUnits: (kpc, 10^10 Msun, Gyr, km/s)

Notes
-----
Astropy's type stubs are incomplete, so we ignore type errors here.
"""

from typing import ClassVar

from astropy import constants, units as u
import numpy as np

_G_const = constants.G
_c_const = constants.c
_M10 = 1e10 * u.Msun


class StellarUnits:
    """Physical constants in stellar units (pc, Msun, Myr, km/s).

    All values are floats (dimensionless) in the specified unit system.

    Attributes
    ----------
    G : float
        Gravitational constant in pc³ / (Msun Myr²).
    G_kms : float
        Gravitational constant in pc (km/s)² / Msun.
    c : float
        Speed of light in pc / Myr.
    H0 : float
        Hubble constant in 1/Myr.
    h : float
        Dimensionless Hubble parameter (H0 / 100 km/s/Mpc).
    rho_crit : float
        Critical density in Msun / pc³.

    Examples
    --------
    >>> units = StellarUnits()
    >>> units.G
    0.004498502151469553
    >>> units = StellarUnits(h=0.674)
    >>> units.h
    0.674
    """

    __slots__ = ("h", "H0", "rho_crit")

    # Class constants (computed once at import)
    G: ClassVar[float] = float(_G_const.to(u.pc**3 / u.Msun / u.Myr**2).value)
    G_kms: ClassVar[float] = float(_G_const.to(u.pc * (u.km / u.s) ** 2 / u.Msun).value)
    c: ClassVar[float] = float(_c_const.to(u.pc / u.Myr).value)

    h: float
    H0: float
    rho_crit: float

    def __init__(self, h: float = 0.7) -> None:
        """Create units with specified Hubble parameter.

        Parameters
        ----------
        h : float, optional
            Dimensionless Hubble parameter (H0 / 100 km/s/Mpc). Default is 0.7.
        """
        H0_quantity = h * 100.0 * u.km / u.s / u.Mpc
        object.__setattr__(self, "h", h)
        object.__setattr__(self, "H0", float(H0_quantity.to(1 / u.Myr).value))
        object.__setattr__(
            self,
            "rho_crit",
            float((3 * H0_quantity**2 / (8 * np.pi * _G_const)).to(u.Msun / u.pc**3).value),
        )

    def __setattr__(self, name: str, value: object) -> None:
        """Prevent modification of attributes."""
        raise AttributeError("StellarUnits instances are immutable")

    def __repr__(self) -> str:
        """Return string representation."""
        return f"StellarUnits(h={self.h})"


class GalacticUnits:
    """Physical constants in galactic units (kpc, 10^10 Msun, Gyr, km/s).

    All values are floats (dimensionless) in the specified unit system.

    Attributes
    ----------
    G : float
        Gravitational constant in kpc³ / (10^10 Msun Gyr²).
    G_kms : float
        Gravitational constant in kpc (km/s)² / (10^10 Msun).
    c : float
        Speed of light in kpc / Gyr.
    H0 : float
        Hubble constant in 1/Gyr.
    h : float
        Dimensionless Hubble parameter (H0 / 100 km/s/Mpc).
    rho_crit : float
        Critical density in 10^10 Msun / kpc³.

    Examples
    --------
    >>> units = GalacticUnits()
    >>> units.G
    44985.02151469553
    >>> units = GalacticUnits(h=0.674)
    >>> units.h
    0.674
    """

    __slots__ = ("h", "H0", "rho_crit")

    # Class constants (computed once at import)
    G: ClassVar[float] = float(_G_const.to(u.kpc**3 / _M10 / u.Gyr**2).value)
    G_kms: ClassVar[float] = float(_G_const.to(u.kpc * (u.km / u.s) ** 2 / _M10).value)
    c: ClassVar[float] = float(_c_const.to(u.kpc / u.Gyr).value)

    h: float
    H0: float
    rho_crit: float

    def __init__(self, h: float = 0.7) -> None:
        """Create units with specified Hubble parameter.

        Parameters
        ----------
        h : float, optional
            Dimensionless Hubble parameter (H0 / 100 km/s/Mpc). Default is 0.7.
        """
        H0_quantity = h * 100.0 * u.km / u.s / u.Mpc
        object.__setattr__(self, "h", h)
        object.__setattr__(self, "H0", float(H0_quantity.to(1 / u.Gyr).value))
        object.__setattr__(
            self,
            "rho_crit",
            float((3 * H0_quantity**2 / (8 * np.pi * _G_const)).to(_M10 / u.kpc**3).value),
        )

    def __setattr__(self, name: str, value: object) -> None:
        """Prevent modification of attributes."""
        raise AttributeError("GalacticUnits instances are immutable")

    def __repr__(self) -> str:
        """Return string representation."""
        return f"GalacticUnits(h={self.h})"
