# pylint: disable=no-member
"""Unit tests for physical constants and units module."""

from astropy import constants, units as u
import numpy as np
import pytest

from galactic_dynamics_bovy.utils.units import GalacticUnits, StellarUnits

# ── Reference values from astropy (used for cross-checks) ────────────────────

_G = constants.G
_c = constants.c
_M10 = 1e10 * u.Msun


class TestStellarUnitsClassConstants:
    """Tests for StellarUnits class-level constants (G, G_kms, c)."""

    def test_G_is_positive(self) -> None:
        """Gravitational constant must be positive."""
        assert StellarUnits.G > 0

    def test_G_matches_astropy(self) -> None:
        """G should match astropy conversion to pc^3/(Msun Myr^2)."""
        expected = float(_G.to(u.pc**3 / u.Msun / u.Myr**2).value)
        assert StellarUnits.G == expected

    def test_G_kms_matches_astropy(self) -> None:
        """G_kms should match astropy conversion to pc (km/s)^2/Msun."""
        expected = float(_G.to(u.pc * (u.km / u.s) ** 2 / u.Msun).value)
        assert StellarUnits.G_kms == expected

    def test_c_matches_astropy(self) -> None:
        """Speed of light should match astropy conversion to pc/Myr."""
        expected = float(_c.to(u.pc / u.Myr).value)
        assert StellarUnits.c == expected

    def test_c_order_of_magnitude(self) -> None:
        """Speed of light in pc/Myr should be ~306,600."""
        assert 3.0e5 < StellarUnits.c < 3.1e5

    def test_class_constants_same_on_instance(self) -> None:
        """Class constants should be accessible from instances."""
        s = StellarUnits()
        assert s.G == StellarUnits.G
        assert s.G_kms == StellarUnits.G_kms
        assert s.c == StellarUnits.c


class TestStellarUnitsInstance:
    """Tests for StellarUnits instance attributes (h, H0, rho_crit)."""

    def test_default_h(self) -> None:
        """Default Hubble parameter should be 0.7."""
        s = StellarUnits()
        assert s.h == 0.7

    def test_custom_h(self) -> None:
        """Custom Hubble parameter should be stored."""
        s = StellarUnits(h=0.674)
        assert s.h == 0.674

    def test_H0_matches_astropy(self) -> None:
        """H0 should match astropy conversion to 1/Myr for default h=0.7."""
        H0_quantity = 0.7 * 100.0 * u.km / u.s / u.Mpc
        expected = float(H0_quantity.to(1 / u.Myr).value)
        s = StellarUnits()
        assert s.H0 == expected

    def test_H0_scales_with_h(self) -> None:
        """H0 should scale linearly with h."""
        s1 = StellarUnits(h=0.7)
        s2 = StellarUnits(h=1.0)
        assert np.isclose(s2.H0 / s1.H0, 1.0 / 0.7)

    def test_rho_crit_is_positive(self) -> None:
        """Critical density must be positive."""
        s = StellarUnits()
        assert s.rho_crit > 0

    def test_rho_crit_scales_with_h_squared(self) -> None:
        """Critical density should scale as h^2."""
        s1 = StellarUnits(h=0.7)
        s2 = StellarUnits(h=1.0)
        assert np.isclose(s2.rho_crit / s1.rho_crit, (1.0 / 0.7) ** 2)


class TestStellarUnitsImmutability:
    """Tests for StellarUnits immutability."""

    def test_cannot_set_h(self) -> None:
        """Setting h should raise AttributeError."""
        s = StellarUnits()
        with pytest.raises(AttributeError, match="immutable"):
            s.h = 0.5  # type: ignore[misc]

    def test_cannot_set_H0(self) -> None:
        """Setting H0 should raise AttributeError."""
        s = StellarUnits()
        with pytest.raises(AttributeError, match="immutable"):
            s.H0 = 1.0  # type: ignore[misc]

    def test_cannot_set_new_attribute(self) -> None:
        """Setting a new attribute should raise AttributeError."""
        s = StellarUnits()
        with pytest.raises(AttributeError):
            s.foo = 42  # type: ignore[attr-defined]

    def test_repr(self) -> None:
        """Repr should show the h parameter."""
        s = StellarUnits(h=0.674)
        assert repr(s) == "StellarUnits(h=0.674)"


# ── GalacticUnits ─────────────────────────────────────────────────────────────


class TestGalacticUnitsClassConstants:
    """Tests for GalacticUnits class-level constants (G, G_kms, c)."""

    def test_G_is_positive(self) -> None:
        """Gravitational constant must be positive."""
        assert GalacticUnits.G > 0

    def test_G_matches_astropy(self) -> None:
        """G should match astropy conversion to kpc^3/(10^10 Msun Gyr^2)."""
        expected = float(_G.to(u.kpc**3 / _M10 / u.Gyr**2).value)
        assert GalacticUnits.G == expected

    def test_G_kms_matches_astropy(self) -> None:
        """G_kms should match astropy conversion to kpc (km/s)^2/(10^10 Msun)."""
        expected = float(_G.to(u.kpc * (u.km / u.s) ** 2 / _M10).value)
        assert GalacticUnits.G_kms == expected

    def test_c_matches_astropy(self) -> None:
        """Speed of light should match astropy conversion to kpc/Gyr."""
        expected = float(_c.to(u.kpc / u.Gyr).value)
        assert GalacticUnits.c == expected

    def test_c_order_of_magnitude(self) -> None:
        """Speed of light in kpc/Gyr should be ~306,600."""
        assert 3.0e5 < GalacticUnits.c < 3.1e5

    def test_class_constants_same_on_instance(self) -> None:
        """Class constants should be accessible from instances."""
        g = GalacticUnits()
        assert g.G == GalacticUnits.G
        assert g.G_kms == GalacticUnits.G_kms
        assert g.c == GalacticUnits.c


class TestGalacticUnitsInstance:
    """Tests for GalacticUnits instance attributes (h, H0, rho_crit)."""

    def test_default_h(self) -> None:
        """Default Hubble parameter should be 0.7."""
        g = GalacticUnits()
        assert g.h == 0.7

    def test_custom_h(self) -> None:
        """Custom Hubble parameter should be stored."""
        g = GalacticUnits(h=0.674)
        assert g.h == 0.674

    def test_H0_matches_astropy(self) -> None:
        """H0 should match astropy conversion to 1/Gyr for default h=0.7."""
        H0_quantity = 0.7 * 100.0 * u.km / u.s / u.Mpc
        expected = float(H0_quantity.to(1 / u.Gyr).value)
        g = GalacticUnits()
        assert g.H0 == expected

    def test_H0_scales_with_h(self) -> None:
        """H0 should scale linearly with h."""
        g1 = GalacticUnits(h=0.7)
        g2 = GalacticUnits(h=1.0)
        assert np.isclose(g2.H0 / g1.H0, 1.0 / 0.7)

    def test_rho_crit_is_positive(self) -> None:
        """Critical density must be positive."""
        g = GalacticUnits()
        assert g.rho_crit > 0

    def test_rho_crit_scales_with_h_squared(self) -> None:
        """Critical density should scale as h^2."""
        g1 = GalacticUnits(h=0.7)
        g2 = GalacticUnits(h=1.0)
        assert np.isclose(g2.rho_crit / g1.rho_crit, (1.0 / 0.7) ** 2)


class TestGalacticUnitsImmutability:
    """Tests for GalacticUnits immutability."""

    def test_cannot_set_h(self) -> None:
        """Setting h should raise AttributeError."""
        g = GalacticUnits()
        with pytest.raises(AttributeError, match="immutable"):
            g.h = 0.5  # type: ignore[misc]

    def test_cannot_set_rho_crit(self) -> None:
        """Setting rho_crit should raise AttributeError."""
        g = GalacticUnits()
        with pytest.raises(AttributeError, match="immutable"):
            g.rho_crit = 1.0  # type: ignore[misc]

    def test_cannot_set_new_attribute(self) -> None:
        """Setting a new attribute should raise AttributeError."""
        g = GalacticUnits()
        with pytest.raises(AttributeError):
            g.foo = 42  # type: ignore[attr-defined]

    def test_repr(self) -> None:
        """Repr should show the h parameter."""
        g = GalacticUnits(h=0.674)
        assert repr(g) == "GalacticUnits(h=0.674)"


# ── Cross-checks between unit systems ────────────────────────────────────────


class TestUnitSystemConsistency:
    """Cross-checks between StellarUnits and GalacticUnits."""

    def test_G_ratio(self) -> None:
        """G_galactic / G_stellar = (pc/kpc)^3 x (M10/Msun) x (Gyr/Myr)^2 = 1e7."""
        ratio = GalacticUnits.G / StellarUnits.G
        # pc^3->kpc^3: (1e-3)^3=1e-9, Msun->M10: 1e10, Myr^2->Gyr^2: (1e3)^2=1e6
        expected = 1e-9 * 1e10 * 1e6
        assert np.isclose(ratio, expected)

    def test_c_is_same_in_both_systems(self) -> None:
        """c in kpc/Gyr should equal c in pc/Myr (both ~306)."""
        # kpc/Gyr = 1000 pc / 1000 Myr = pc/Myr
        assert np.isclose(GalacticUnits.c, StellarUnits.c)

    def test_rho_crit_ratio(self) -> None:
        """rho_crit ratio should reflect mass/length^3 unit change."""
        s = StellarUnits()
        g = GalacticUnits()
        # Msun/pc^3 -> (10^10 Msun)/kpc^3 requires factor 1e-10 * (1e3)^3 = 1e-1
        ratio = s.rho_crit / g.rho_crit
        expected = 1e10 / (1e3) ** 3  # Msun/(10^10 Msun) * kpc^3/pc^3
        assert np.isclose(ratio, expected)
