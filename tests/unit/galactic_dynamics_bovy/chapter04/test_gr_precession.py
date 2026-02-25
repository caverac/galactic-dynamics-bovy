"""Unit tests for GR apsidal precession module."""

from pathlib import Path
from unittest.mock import MagicMock, patch

import numpy as np

from galactic_dynamics_bovy.chapter04.gr_precession import (
    _integrand,
    _setup_axes,
    _solve_EL,
    compute_delta_psi,
    compute_precession_curve,
    plot_gr_precession,
)
from galactic_dynamics_bovy.utils.units import SolarUnits

_GM_SUN = SolarUnits.GM
_C = SolarUnits.c
_AU = SolarUnits.AU


class TestSolveEL:
    """Tests for _solve_EL."""

    def test_turning_points_satisfy_energy(self) -> None:
        """E should equal Phi_eff at both r_p and r_a."""
        r_p = 0.3 * _AU
        r_a = 0.5 * _AU
        E, L2 = _solve_EL(r_p, r_a, _GM_SUN, _C)
        c2 = _C**2
        phi_p = -_GM_SUN / r_p + L2 / (2.0 * r_p**2) - _GM_SUN * L2 / (c2 * r_p**3)
        phi_a = -_GM_SUN / r_a + L2 / (2.0 * r_a**2) - _GM_SUN * L2 / (c2 * r_a**3)
        assert np.isclose(E, phi_p, rtol=1e-12)
        assert np.isclose(E, phi_a, rtol=1e-12)

    def test_negative_energy(self) -> None:
        """Bound orbit should have negative energy."""
        r_p = 0.3 * _AU
        r_a = 0.5 * _AU
        E, _ = _solve_EL(r_p, r_a, _GM_SUN, _C)
        assert E < 0

    def test_positive_L2(self) -> None:
        """Squared angular momentum should be positive."""
        r_p = 0.3 * _AU
        r_a = 0.5 * _AU
        _, L2 = _solve_EL(r_p, r_a, _GM_SUN, _C)
        assert L2 > 0

    def test_approaches_keplerian(self) -> None:
        """For large c, L² should approach the Keplerian value GM*r_p*r_a/(r_p+r_a)/0.5."""
        r_p = 0.3 * _AU
        r_a = 0.5 * _AU
        # Keplerian: L² = 2*GM*r_p*r_a / (r_p + r_a)
        L2_kepler = 2.0 * _GM_SUN * r_p * r_a / (r_p + r_a)
        _, L2 = _solve_EL(r_p, r_a, _GM_SUN, 1e20)
        assert np.isclose(L2, L2_kepler, rtol=1e-8)


class TestComputeDeltaPsi:
    """Tests for compute_delta_psi."""

    def test_mercury_precession_analytical(self) -> None:
        """GR precession for Mercury should match the analytical prediction.

        Analytical: delta_phi = 6*pi*GM / (c^2 * r_p * (1+e)).
        """
        e = 0.205630
        a = 0.387098 * _AU
        r_p = a * (1.0 - e)
        analytical = 6.0 * np.pi * _GM_SUN / (_C**2 * r_p * (1.0 + e))
        numerical = compute_delta_psi(r_p, e) - 2.0 * np.pi
        assert np.isclose(numerical, analytical, rtol=1e-3)

    def test_mercury_precession_43_arcsec_century(self) -> None:
        """GR precession for Mercury should be ~43 arcsec/century."""
        e = 0.205630
        a = 0.387098 * _AU
        r_p = a * (1.0 - e)
        precession_rad = compute_delta_psi(r_p, e) - 2.0 * np.pi
        precession_arcsec = np.degrees(precession_rad) * 3600.0
        T_mercury = 2.0 * np.pi * np.sqrt(a**3 / _GM_SUN)  # seconds
        T_century = 100.0 * 365.25 * 86400.0  # seconds
        arcsec_per_century = precession_arcsec * T_century / T_mercury
        assert np.isclose(arcsec_per_century, 43.0, atol=0.2)

    def test_near_circular(self) -> None:
        """Near-circular orbit should have precession close to 6*pi*GM/(c^2*r_p*(1+e))."""
        e = 0.05
        r_p = 0.3 * _AU
        analytical = 6.0 * np.pi * _GM_SUN / (_C**2 * r_p * (1.0 + e))
        numerical = compute_delta_psi(r_p, e) - 2.0 * np.pi
        assert np.isclose(numerical, analytical, rtol=1e-3)

    def test_positive_precession(self) -> None:
        """GR precession should be prograde (positive)."""
        r_p = 0.3 * _AU
        e = 0.3
        delta = compute_delta_psi(r_p, e) - 2.0 * np.pi
        assert delta > 0

    def test_scales_with_inverse_rp(self) -> None:
        """At fixed e, precession should scale roughly as 1/r_p."""
        e = 0.2
        r1 = 0.2 * _AU
        r2 = 0.4 * _AU
        d1 = compute_delta_psi(r1, e) - 2.0 * np.pi
        d2 = compute_delta_psi(r2, e) - 2.0 * np.pi
        ratio = d1 / d2
        assert np.isclose(ratio, 2.0, rtol=0.01)


class TestComputePrecessionCurve:
    """Tests for compute_precession_curve."""

    def test_output_shape(self) -> None:
        """Output should match input array length."""
        a = np.array([0.2, 0.4, 0.6]) * _AU
        result = compute_precession_curve(a, 0.2)
        assert result.shape == (3,)

    def test_all_positive(self) -> None:
        """All precession values should be positive for GR correction."""
        a = np.linspace(0.2, 1.0, 5) * _AU
        result = compute_precession_curve(a, 0.2)
        assert np.all(result > 0)

    def test_decreasing_with_a(self) -> None:
        """Precession should decrease as semi-major axis increases."""
        a = np.linspace(0.2, 1.0, 5) * _AU
        result = compute_precession_curve(a, 0.2)
        assert np.all(np.diff(result) < 0)


class TestIntegrand:
    """Tests for _integrand helper."""

    def test_returns_zero_when_kinetic_negative(self) -> None:
        """Guard should return 0 when E < Phi_eff."""
        # k0 extremely negative guarantees kinetic < 0
        result = _integrand(0.0, 2.0, 1.0, 1.0, (-1e30, 0.0, 0.0, 0.0))
        assert result == 0.0

    def test_positive_for_valid_orbit(self) -> None:
        """Should return a positive value at t=0 for a valid orbit."""
        r_p = 0.3 * _AU
        r_a = 0.5 * _AU
        E, L2 = _solve_EL(r_p, r_a, _GM_SUN, _C)
        S = r_p + r_a
        D = r_a - r_p
        LD = np.sqrt(L2) * D
        c2 = _C**2
        coeffs = (2.0 * E, 2.0 * _GM_SUN, -L2, 2.0 * _GM_SUN * L2 / c2)
        val = _integrand(0.0, S, D, LD, coeffs)
        assert val > 0


class TestSetupAxes:
    """Tests for _setup_axes helper."""

    def test_applies_formatting(self) -> None:
        """Should call tick_params and grid on the axes."""
        mock_ax = MagicMock()
        _setup_axes(mock_ax)
        assert mock_ax.tick_params.called
        mock_ax.grid.assert_called_once()


class TestPlotGrPrecession:
    """Tests for plot_gr_precession function."""

    @patch("galactic_dynamics_bovy.chapter04.gr_precession.compute_precession_curve")
    @patch("galactic_dynamics_bovy.chapter04.gr_precession.plt")
    def test_creates_figure(self, mock_plt: MagicMock, mock_curve: MagicMock) -> None:
        """Verify a figure is created."""
        mock_fig = MagicMock()
        mock_ax = MagicMock()
        mock_plt.subplots.return_value = (mock_fig, mock_ax)
        mock_curve.return_value = np.ones(200) * 1e-7
        plot_gr_precession()
        mock_plt.subplots.assert_called_once()

    @patch("galactic_dynamics_bovy.chapter04.gr_precession.compute_precession_curve")
    @patch("galactic_dynamics_bovy.chapter04.gr_precession.plt")
    def test_shows_plot_when_no_path(self, mock_plt: MagicMock, mock_curve: MagicMock) -> None:
        """Verify plt.show() is called when no path provided."""
        mock_plt.subplots.return_value = (MagicMock(), MagicMock())
        mock_curve.return_value = np.ones(200) * 1e-7
        plot_gr_precession(path=None)
        mock_plt.show.assert_called_once()

    @patch("galactic_dynamics_bovy.chapter04.gr_precession.compute_precession_curve")
    @patch("galactic_dynamics_bovy.chapter04.gr_precession.plt")
    def test_saves_figure_when_path_provided(self, mock_plt: MagicMock, mock_curve: MagicMock) -> None:
        """Verify figure is saved when path provided."""
        mock_fig = MagicMock()
        mock_ax = MagicMock()
        mock_plt.subplots.return_value = (mock_fig, mock_ax)
        mock_curve.return_value = np.ones(200) * 1e-7
        test_path = Path("/tmp/test_gr_precession.png")
        plot_gr_precession(path=test_path)
        mock_fig.savefig.assert_called_once()
        mock_plt.close.assert_called_once_with(mock_fig)

    @patch("galactic_dynamics_bovy.chapter04.gr_precession.compute_precession_curve")
    @patch("galactic_dynamics_bovy.chapter04.gr_precession.plt")
    def test_plots_three_curves(self, mock_plt: MagicMock, mock_curve: MagicMock) -> None:
        """Should call compute_precession_curve three times (one per eccentricity)."""
        mock_plt.subplots.return_value = (MagicMock(), MagicMock())
        mock_curve.return_value = np.ones(200) * 1e-7
        plot_gr_precession()
        assert mock_curve.call_count == 3

    @patch("galactic_dynamics_bovy.chapter04.gr_precession.compute_precession_curve")
    @patch("galactic_dynamics_bovy.chapter04.gr_precession.plt")
    def test_sets_labels(self, mock_plt: MagicMock, mock_curve: MagicMock) -> None:
        """Verify axis labels are set."""
        mock_ax = MagicMock()
        mock_plt.subplots.return_value = (MagicMock(), mock_ax)
        mock_curve.return_value = np.ones(200) * 1e-7
        plot_gr_precession()
        mock_ax.set_xlabel.assert_called_once()
        mock_ax.set_ylabel.assert_called_once()
        mock_ax.legend.assert_called_once()
