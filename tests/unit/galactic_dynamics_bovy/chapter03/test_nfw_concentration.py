"""Unit tests for NFW concentration module."""

from unittest.mock import patch

import numpy as np
import pytest

from galactic_dynamics_bovy.chapter03.nfw_concentration import (
    enclosed_mass_residual,
    m200_from_concentration,
    nfw_f,
    nfw_potential,
    nfw_potential_exterior,
    nfw_potential_interior,
    nfw_vesc_to_r200,
    r200_from_m200,
    solve_concentration,
    solve_concentration_from_vesc,
    SolvedNFWModel,
)
from galactic_dynamics_bovy.utils.units import GalacticUnits


class TestNfwF:
    """Tests for nfw_f function."""

    def test_zero_at_zero(self) -> None:
        """f(0) should be 0."""
        assert np.isclose(nfw_f(0.0), 0.0, atol=1e-10)

    def test_positive_for_positive_x(self) -> None:
        """f(x) should be positive for x > 0."""
        for x in [0.1, 1.0, 10.0, 100.0]:
            assert nfw_f(x) > 0

    def test_monotonically_increasing(self) -> None:
        """f(x) should be monotonically increasing."""
        x_values = [0.1, 1.0, 5.0, 10.0, 20.0]
        f_values = [nfw_f(x) for x in x_values]
        assert all(f_values[i] < f_values[i + 1] for i in range(len(f_values) - 1))

    def test_known_value(self) -> None:
        """Test against known value at x=1."""
        # f(1) = ln(2) - 1/2 ~ 0.193
        expected = np.log(2) - 0.5
        assert np.isclose(nfw_f(1.0), expected)


class TestM200FromConcentration:
    """Tests for m200_from_concentration function."""

    def test_returns_positive(self) -> None:
        """M_200 should be positive."""
        for c in [5.0, 10.0, 15.0]:
            assert m200_from_concentration(c) > 0

    def test_decreases_with_concentration(self) -> None:
        """Higher concentration should give lower M_200 (D&M relation)."""
        m1 = m200_from_concentration(5.0)
        m2 = m200_from_concentration(10.0)
        m3 = m200_from_concentration(15.0)
        assert m1 > m2 > m3

    def test_known_value(self) -> None:
        """Test against known value from D&M relation."""
        # For c = 10^0.905 ~ 8.03, M_200 = 10^12 h^-1 M_sun
        c_ref = 10**0.905
        m200 = m200_from_concentration(c_ref, h=0.7)
        # Expected: 10^12 / 0.7 M_sun = 100 / 0.7 x 10^10 M_sun ~ 143
        expected = 100.0 / 0.7
        assert np.isclose(m200, expected, rtol=0.01)

    def test_h_dependence(self) -> None:
        """M_200 should scale with 1/h."""
        m1 = m200_from_concentration(10.0, h=0.7)
        m2 = m200_from_concentration(10.0, h=1.0)
        assert np.isclose(m1 / m2, 1.0 / 0.7, rtol=0.01)


class TestR200FromM200:
    """Tests for r200_from_m200 function."""

    def test_returns_positive(self) -> None:
        """R_200 should be positive."""
        units = GalacticUnits()
        assert r200_from_m200(100.0, units.rho_crit) > 0

    def test_scales_as_cube_root(self) -> None:
        """R_200 should scale as M^(1/3)."""
        units = GalacticUnits()
        r1 = r200_from_m200(100.0, units.rho_crit)
        r2 = r200_from_m200(800.0, units.rho_crit)  # 8x mass
        assert np.isclose(r2 / r1, 2.0, rtol=0.01)  # 2x radius

    def test_typical_value(self) -> None:
        """Test for MW-like halo."""
        units = GalacticUnits(h=0.7)
        # M_200 = 10^12 M_sun = 100 x 10^10 M_sun
        r200 = r200_from_m200(100.0, units.rho_crit)
        # R_200 should be ~200 kpc for MW-like halo
        assert 150 < r200 < 250


class TestEnclosedMassResidual:
    """Tests for enclosed_mass_residual function."""

    def test_zero_at_solution(self) -> None:
        """Residual should be zero at the solution."""
        result = solve_concentration(h=0.7)
        units = GalacticUnits(h=0.7)

        residual = enclosed_mass_residual(result.c, r_inner=8.0, m_inner=9.0, rho_crit=units.rho_crit, h=0.7)
        assert np.isclose(residual, 0.0, atol=1e-6)

    def test_sign_change(self) -> None:
        """Residual should change sign across the solution."""
        units = GalacticUnits(h=0.7)

        res_low = enclosed_mass_residual(5.0, r_inner=8.0, m_inner=9.0, rho_crit=units.rho_crit, h=0.7)
        res_high = enclosed_mass_residual(6.0, r_inner=8.0, m_inner=9.0, rho_crit=units.rho_crit, h=0.7)
        # Should have opposite signs (or one is zero)
        assert res_low * res_high < 0


class TestSolveConcentration:
    """Tests for solve_concentration function."""

    def test_returns_dataclass(self) -> None:
        """Should return SolvedNFWModel dataclass."""
        result = solve_concentration()
        assert isinstance(result, SolvedNFWModel)
        assert hasattr(result, "c")
        assert hasattr(result, "m200")
        assert hasattr(result, "r200")
        assert hasattr(result, "a")

    def test_all_positive(self) -> None:
        """All values should be positive."""
        result = solve_concentration()
        assert result.c > 0
        assert result.m200 > 0
        assert result.r200 > 0
        assert result.a > 0

    def test_consistency(self) -> None:
        """a = R_200 / c should hold."""
        result = solve_concentration()
        assert np.isclose(result.a, result.r200 / result.c)

    def test_satisfies_constraint(self) -> None:
        """Solution should satisfy the enclosed mass constraint."""
        r_inner = 8.0
        m_inner = 9.0
        result = solve_concentration(r_inner=r_inner, m_inner=m_inner)

        # Verify M(<r_inner) = m_inner
        x = r_inner * result.c / result.r200
        m_enclosed = result.m200 * nfw_f(x) / nfw_f(result.c)
        assert np.isclose(m_enclosed, m_inner, rtol=1e-6)

    def test_different_constraints(self) -> None:
        """Different m_inner should give different solutions."""
        result1 = solve_concentration(m_inner=9.0)
        result2 = solve_concentration(m_inner=5.0)
        assert result1.c != result2.c
        assert result1.m200 != result2.m200

    def test_different_h(self) -> None:
        """Different h should give different solutions."""
        result1 = solve_concentration(h=0.7)
        result2 = solve_concentration(h=1.0)
        assert result1.c != result2.c

    def test_no_solution_raises(self) -> None:
        """Should raise if no solution in bracket."""
        with pytest.raises(ValueError):
            # Very small enclosed mass - no solution in bracket
            solve_concentration(m_inner=0.1, c_bracket=(5.0, 6.0))


class TestNfwPotentialInterior:
    """Tests for nfw_potential_interior function."""

    def test_negative(self) -> None:
        """Interior potential should be negative."""
        units = GalacticUnits()
        phi = nfw_potential_interior(r=10.0, m200=100.0, c=10.0, r200=200.0, G=units.G_kms)
        assert phi < 0

    def test_equals_point_mass(self) -> None:
        """Interior potential should equal -G M(<r) / r."""
        units = GalacticUnits()
        r, m200, c, r200 = 10.0, 100.0, 10.0, 200.0
        a = r200 / c
        fc = nfw_f(c)
        m_enclosed = m200 * nfw_f(r / a) / fc
        expected = -units.G_kms * m_enclosed / r
        actual = nfw_potential_interior(r, m200, c, r200, units.G_kms)
        assert np.isclose(actual, expected)


class TestNfwPotentialExterior:
    """Tests for nfw_potential_exterior function."""

    def test_negative(self) -> None:
        """Exterior potential should be negative (inside the shell)."""
        units = GalacticUnits()
        phi = nfw_potential_exterior(r=10.0, m200=100.0, c=10.0, r200=200.0, G=units.G_kms)
        assert phi < 0

    def test_zero_at_r200(self) -> None:
        """Exterior potential should be zero at R_200 (no mass outside)."""
        units = GalacticUnits()
        phi = nfw_potential_exterior(r=200.0, m200=100.0, c=10.0, r200=200.0, G=units.G_kms)
        assert np.isclose(phi, 0.0, atol=1e-10)

    def test_increases_with_radius(self) -> None:
        """Exterior potential should become less negative with radius."""
        units = GalacticUnits()
        phi1 = nfw_potential_exterior(r=10.0, m200=100.0, c=10.0, r200=200.0, G=units.G_kms)
        phi2 = nfw_potential_exterior(r=50.0, m200=100.0, c=10.0, r200=200.0, G=units.G_kms)
        phi3 = nfw_potential_exterior(r=100.0, m200=100.0, c=10.0, r200=200.0, G=units.G_kms)
        assert phi1 < phi2 < phi3


class TestNfwPotential:
    """Tests for nfw_potential function."""

    def test_negative(self) -> None:
        """Potential should be negative."""
        units = GalacticUnits()
        phi = nfw_potential(r=10.0, m200=100.0, c=10.0, r200=200.0, G=units.G_kms)
        assert phi < 0

    def test_sum_of_interior_and_exterior(self) -> None:
        """Total potential should equal interior + exterior."""
        units = GalacticUnits()
        r, m200, c, r200 = 10.0, 100.0, 10.0, 200.0
        phi_int = nfw_potential_interior(r, m200, c, r200, units.G_kms)
        phi_ext = nfw_potential_exterior(r, m200, c, r200, units.G_kms)
        phi_total = nfw_potential(r, m200, c, r200, units.G_kms)
        assert np.isclose(phi_total, phi_int + phi_ext)

    def test_increases_with_radius(self) -> None:
        """Potential should become less negative with increasing radius."""
        units = GalacticUnits()
        phi1 = nfw_potential(r=10.0, m200=100.0, c=10.0, r200=200.0, G=units.G_kms)
        phi2 = nfw_potential(r=50.0, m200=100.0, c=10.0, r200=200.0, G=units.G_kms)
        phi3 = nfw_potential(r=100.0, m200=100.0, c=10.0, r200=200.0, G=units.G_kms)
        assert phi1 < phi2 < phi3

    def test_scales_with_mass(self) -> None:
        """Potential should scale linearly with mass."""
        units = GalacticUnits()
        phi1 = nfw_potential(r=10.0, m200=100.0, c=10.0, r200=200.0, G=units.G_kms)
        phi2 = nfw_potential(r=10.0, m200=200.0, c=10.0, r200=200.0, G=units.G_kms)
        assert np.isclose(phi2, 2 * phi1)


class TestNfwVescToR200:
    """Tests for nfw_vesc_to_r200 function."""

    def test_positive(self) -> None:
        """Escape velocity should be positive."""
        units = GalacticUnits()
        v = nfw_vesc_to_r200(r=10.0, m200=100.0, c=10.0, r200=200.0, G=units.G_kms)
        assert v > 0

    def test_decreases_with_radius(self) -> None:
        """Escape velocity should decrease with radius (to fixed R_200)."""
        units = GalacticUnits()
        v1 = nfw_vesc_to_r200(r=10.0, m200=100.0, c=10.0, r200=200.0, G=units.G_kms)
        v2 = nfw_vesc_to_r200(r=50.0, m200=100.0, c=10.0, r200=200.0, G=units.G_kms)
        v3 = nfw_vesc_to_r200(r=100.0, m200=100.0, c=10.0, r200=200.0, G=units.G_kms)
        assert v1 > v2 > v3

    def test_zero_at_r200(self) -> None:
        """Escape velocity should be zero at R_200 (escaping to R_200)."""
        units = GalacticUnits()
        v = nfw_vesc_to_r200(r=200.0, m200=100.0, c=10.0, r200=200.0, G=units.G_kms)
        assert np.isclose(v, 0.0, atol=1e-10)


class TestSolveConcentrationFromVesc:
    """Tests for solve_concentration_from_vesc function."""

    def test_returns_dataclass(self) -> None:
        """Should return SolvedNFWModel dataclass."""
        result = solve_concentration_from_vesc()
        assert isinstance(result, SolvedNFWModel)
        assert hasattr(result, "c")
        assert hasattr(result, "m200")
        assert hasattr(result, "r200")
        assert hasattr(result, "a")
        assert hasattr(result, "v_esc_predicted")

    def test_all_positive(self) -> None:
        """All values should be positive."""
        result = solve_concentration_from_vesc()
        assert result.c > 0
        assert result.m200 > 0
        assert result.r200 > 0
        assert result.a > 0
        assert result.v_esc_predicted > 0

    def test_consistency(self) -> None:
        """a = R_200 / c should hold."""
        result = solve_concentration_from_vesc()
        assert np.isclose(result.a, result.r200 / result.c)

    def test_satisfies_mass_constraint(self) -> None:
        """Solution should satisfy the enclosed mass constraint."""
        r_inner = 8.0
        m_inner = 9.0
        result = solve_concentration_from_vesc(r_inner=r_inner, m_inner=m_inner)

        # Verify M(<r_inner) = m_inner
        x = r_inner / result.a
        m_enclosed = result.m200 * nfw_f(x) / nfw_f(result.c)
        assert np.isclose(m_enclosed, m_inner, rtol=1e-4)

    def test_satisfies_vesc_constraint(self) -> None:
        """Solution should satisfy the escape velocity constraint."""
        v_esc = 550.0
        result = solve_concentration_from_vesc(v_esc=v_esc)
        assert np.isclose(result.v_esc_predicted, v_esc, rtol=1e-4)

    def test_reasonable_mw_mass(self) -> None:
        """Solution should give reasonable MW mass (~10^12 M_sun)."""
        result = solve_concentration_from_vesc(m_inner=9.0, v_esc=550.0)
        # M_200 should be roughly 10^12 M_sun = 100 x 10^10 M_sun
        assert 50 < result.m200 < 300  # Allow factor of ~3 uncertainty

    def test_different_vesc_gives_different_solution(self) -> None:
        """Different v_esc should give different solutions."""
        result1 = solve_concentration_from_vesc(v_esc=500.0)
        result2 = solve_concentration_from_vesc(v_esc=600.0)
        assert result1.c != result2.c
        assert result1.m200 != result2.m200

    def test_convergence_failure_raises(self) -> None:
        """Should raise RuntimeError if solver fails to converge."""
        # Mock fsolve to return a failed convergence status (ier != 1)
        mock_result = (
            np.array([10.0, 100.0]),  # solution
            {},  # info dict
            0,  # ier = 0 means failed
            "Mocked convergence failure",  # message
        )
        with patch("galactic_dynamics_bovy.chapter03.nfw_concentration.optimize.fsolve") as mock_fsolve:
            mock_fsolve.return_value = mock_result
            with pytest.raises(RuntimeError, match="Failed to converge"):
                solve_concentration_from_vesc()
