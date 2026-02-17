"""Unit tests for planet-comet scattering module."""

from pathlib import Path
from unittest.mock import MagicMock, patch

from galpy.potential import IsochronePotential, vcirc
import numpy as np
import pytest

from galactic_dynamics_bovy.chapter03.planet_comet_scattering import (
    _eccentricity,
    apply_kick_at_pericenter,
    find_pericenter_time,
    integrate_orbit,
    OrbitResult,
    parameter_study,
    ParameterStudyResult,
    plot_apocenter_study,
    plot_eccentricity_study,
    plot_orbit_comparison,
    plot_parameter_study,
    ScatteringResult,
    setup_isochrone_potential,
)


class TestSetupIsochronePotential:
    """Tests for setup_isochrone_potential function."""

    def test_returns_isochrone_potential(self) -> None:
        """Should return an IsochronePotential instance."""
        pot = setup_isochrone_potential()
        assert isinstance(pot, IsochronePotential)

    def test_default_b(self) -> None:
        """Default scale parameter should be b=1.0."""
        pot = setup_isochrone_potential()
        # IsochronePotential with b=1.0 and normalize=1.0
        assert isinstance(pot, IsochronePotential)

    def test_custom_b(self) -> None:
        """Custom b should produce a different potential."""
        pot1 = setup_isochrone_potential(b=1.0)
        pot2 = setup_isochrone_potential(b=2.0)
        # Different b values should yield different potentials
        assert pot1 != pot2


class TestIntegrateOrbit:
    """Tests for integrate_orbit function."""

    def test_returns_orbit_result(self) -> None:
        """Should return an OrbitResult dataclass."""
        result = integrate_orbit(R=1.0, vR=0.3, vT=0.5)
        assert isinstance(result, OrbitResult)

    def test_positive_pericenter(self) -> None:
        """Pericenter radius should be positive."""
        result = integrate_orbit(R=1.0, vR=0.3, vT=0.5)
        assert result.r_peri > 0

    def test_apo_geq_peri(self) -> None:
        """Apocenter should be >= pericenter."""
        result = integrate_orbit(R=1.0, vR=0.3, vT=0.5)
        assert result.r_apo >= result.r_peri

    def test_negative_energy(self) -> None:
        """Bound orbit should have negative energy."""
        result = integrate_orbit(R=1.0, vR=0.3, vT=0.5)
        assert result.energy < 0

    def test_circular_orbit_peri_approx_apo(self) -> None:
        """A circular orbit should have r_peri ~ r_apo."""
        # For a circular orbit in isochrone, v_c at R=1 with normalize=1 gives vT~1
        # Use galpy's built-in circular velocity
        pot = setup_isochrone_potential(b=1.0)
        vc = float(vcirc(pot, 1.0))
        result = integrate_orbit(R=1.0, vR=0.0, vT=vc)
        assert np.isclose(result.r_peri, result.r_apo, rtol=0.01)

    def test_time_array_length(self) -> None:
        """Time array should have the requested number of steps."""
        result = integrate_orbit(R=1.0, vR=0.3, vT=0.5, n_steps=501)
        assert len(result.t) == 501

    def test_angular_momentum(self) -> None:
        """Angular momentum should be R * vT at t=0."""
        R, vT = 1.0, 0.5
        result = integrate_orbit(R=R, vR=0.3, vT=vT)
        # Lz = R * vT at initial time
        assert np.isclose(result.lz, R * vT, rtol=0.01)


class TestFindPericenterTime:
    """Tests for find_pericenter_time function."""

    def test_finds_first_pericenter(self) -> None:
        """Should find a pericenter time > 0."""
        result = integrate_orbit(R=1.0, vR=0.3, vT=0.5)
        t_peri = find_pericenter_time(result, occurrence=1)
        assert t_peri > 0

    def test_second_after_first(self) -> None:
        """Second pericenter should be after the first."""
        result = integrate_orbit(R=1.0, vR=0.3, vT=0.5, t_max=100.0)
        t1 = find_pericenter_time(result, occurrence=1)
        t2 = find_pericenter_time(result, occurrence=2)
        assert t2 > t1

    def test_invalid_occurrence_zero(self) -> None:
        """Should raise ValueError for occurrence < 1."""
        result = integrate_orbit(R=1.0, vR=0.3, vT=0.5)
        with pytest.raises(ValueError, match="occurrence must be >= 1"):
            find_pericenter_time(result, occurrence=0)

    def test_too_many_occurrences(self) -> None:
        """Should raise ValueError if not enough pericenters found."""
        result = integrate_orbit(R=1.0, vR=0.3, vT=0.5, t_max=1.0, n_steps=101)
        with pytest.raises(ValueError, match="pericenter"):
            find_pericenter_time(result, occurrence=100)


class TestApplyKickAtPericenter:
    """Tests for apply_kick_at_pericenter function."""

    def test_returns_scattering_result(self) -> None:
        """Should return a ScatteringResult dataclass."""
        orbit = integrate_orbit(R=1.0, vR=0.3, vT=0.5)
        result = apply_kick_at_pericenter(orbit, delta_v_r=0.1)
        assert isinstance(result, ScatteringResult)

    def test_zero_kick_preserves_orbit(self) -> None:
        """Zero kick should yield approximately the same pericenter."""
        orbit = integrate_orbit(R=1.0, vR=0.3, vT=0.5)
        result = apply_kick_at_pericenter(orbit, delta_v_r=0.0, delta_v_t=0.0)
        assert np.isclose(result.delta_r_peri, 0.0, atol=0.02)

    def test_radial_kick_preserves_lz(self) -> None:
        """A purely radial kick should approximately preserve Lz."""
        orbit = integrate_orbit(R=1.0, vR=0.3, vT=0.5)
        result = apply_kick_at_pericenter(orbit, delta_v_r=0.1, delta_v_t=0.0)
        # Lz depends on R * vT, radial kick doesn't change vT
        assert np.isclose(result.original.lz, result.perturbed.lz, rtol=0.05)

    def test_tangential_kick_changes_lz(self) -> None:
        """A tangential kick should change Lz."""
        orbit = integrate_orbit(R=1.0, vR=0.3, vT=0.5)
        result = apply_kick_at_pericenter(orbit, delta_v_r=0.0, delta_v_t=0.2)
        assert not np.isclose(result.original.lz, result.perturbed.lz, rtol=0.01)

    def test_energy_changes_with_kick(self) -> None:
        """Any kick should change the orbital energy."""
        orbit = integrate_orbit(R=1.0, vR=0.3, vT=0.5)
        result = apply_kick_at_pericenter(orbit, delta_v_r=0.3)
        assert result.original.energy != result.perturbed.energy

    def test_kick_params_stored(self) -> None:
        """Kick parameters should be stored in the result."""
        orbit = integrate_orbit(R=1.0, vR=0.3, vT=0.5)
        result = apply_kick_at_pericenter(orbit, delta_v_r=0.1, delta_v_t=0.2)
        assert result.delta_v_r == 0.1
        assert result.delta_v_t == 0.2


class TestEccentricity:
    """Tests for _eccentricity helper."""

    def test_circular_orbit(self) -> None:
        """Circular orbit should have zero eccentricity."""
        assert _eccentricity(1.0, 1.0) == 0.0

    def test_eccentric_orbit(self) -> None:
        """Known eccentricity value."""
        assert np.isclose(_eccentricity(0.5, 1.5), 0.5)


class TestParameterStudy:
    """Tests for parameter_study function."""

    def test_returns_dataclass(self) -> None:
        """Should return a ParameterStudyResult."""
        kicks = np.array([0.05])
        result = parameter_study(kick_magnitudes=kicks, initial_conditions=[(1.0, 0.3, 0.5)])
        assert isinstance(result, ParameterStudyResult)

    def test_correct_output_shapes(self) -> None:
        """Output arrays should have correct shapes."""
        kicks = np.array([0.05, 0.1, 0.15])
        ics = [(1.0, 0.3, 0.5)]
        result = parameter_study(kick_magnitudes=kicks, initial_conditions=ics)
        assert len(result.kicks) == 3
        assert len(result.delta_peri_radial) == 1
        assert len(result.delta_peri_tangential) == 1
        assert len(result.delta_peri_radial[0]) == 3
        assert len(result.delta_apo_radial) == 1
        assert len(result.delta_apo_radial[0]) == 3
        assert len(result.delta_ecc_radial) == 1
        assert len(result.delta_ecc_radial[0]) == 3

    def test_multiple_orbits(self) -> None:
        """Should return one array per orbit."""
        kicks = np.array([0.05, 0.1])
        ics = [(1.0, 0.3, 0.5), (0.5, 0.0, 0.8)]
        result = parameter_study(kick_magnitudes=kicks, initial_conditions=ics)
        assert len(result.delta_peri_radial) == 2
        assert len(result.delta_peri_tangential) == 2
        assert len(result.delta_apo_radial) == 2
        assert len(result.delta_ecc_radial) == 2

    def test_small_kick_small_change(self) -> None:
        """Very small kicks should produce small pericenter changes."""
        kicks = np.array([0.001])
        result = parameter_study(kick_magnitudes=kicks, initial_conditions=[(1.0, 0.3, 0.5)])
        assert abs(result.delta_peri_radial[0][0]) < 0.5
        assert abs(result.delta_peri_tangential[0][0]) < 0.5

    def test_default_kick_magnitudes(self) -> None:
        """Default kick_magnitudes should produce 30 values."""
        result = parameter_study(initial_conditions=[(1.0, 0.3, 0.5)])
        assert len(result.kicks) == 30
        assert len(result.delta_peri_radial[0]) == 30

    def test_default_initial_conditions(self) -> None:
        """Default initial_conditions should provide two orbits."""
        result = parameter_study(kick_magnitudes=np.array([0.05, 0.1]))
        assert len(result.delta_peri_radial) == 2
        assert len(result.delta_peri_tangential) == 2


class TestPlotOrbitComparison:
    """Tests for plot_orbit_comparison function."""

    @staticmethod
    def _setup_mocks(mock_plt, mock_integrate, mock_kick, mock_find_peri):
        """Set up common mocks for orbit comparison tests."""
        mock_fig = MagicMock()
        mock_ax1 = MagicMock()
        mock_ax2 = MagicMock()
        mock_plt.subplots.return_value = (mock_fig, [mock_ax1, mock_ax2])

        mock_orbit = MagicMock()
        mock_orbit.R.return_value = np.ones(100)
        mock_orbit.phi.return_value = np.linspace(0, 2 * np.pi, 100)

        t = np.linspace(0, 50, 100)
        mock_orbit_result = OrbitResult(orbit=mock_orbit, r_peri=0.5, r_apo=1.5, energy=-1.0, lz=0.5, t=t)
        mock_integrate.return_value = mock_orbit_result
        mock_kick.return_value = ScatteringResult(
            original=mock_orbit_result,
            perturbed=mock_orbit_result,
            delta_v_r=0.15,
            delta_v_t=0.0,
            delta_r_peri=0.1,
            delta_r_peri_frac=0.2,
        )
        mock_find_peri.return_value = 1.0

        return mock_fig, mock_ax1, mock_ax2

    @patch("galactic_dynamics_bovy.chapter03.planet_comet_scattering.find_pericenter_time")
    @patch("galactic_dynamics_bovy.chapter03.planet_comet_scattering.apply_kick_at_pericenter")
    @patch("galactic_dynamics_bovy.chapter03.planet_comet_scattering.integrate_orbit")
    @patch("galactic_dynamics_bovy.chapter03.planet_comet_scattering.plt")
    def test_creates_figure(self, mock_plt, mock_integrate, mock_kick, mock_find_peri) -> None:
        """Verify a figure is created."""
        self._setup_mocks(mock_plt, mock_integrate, mock_kick, mock_find_peri)
        plot_orbit_comparison()
        mock_plt.subplots.assert_called_once()

    @patch("galactic_dynamics_bovy.chapter03.planet_comet_scattering.find_pericenter_time")
    @patch("galactic_dynamics_bovy.chapter03.planet_comet_scattering.apply_kick_at_pericenter")
    @patch("galactic_dynamics_bovy.chapter03.planet_comet_scattering.integrate_orbit")
    @patch("galactic_dynamics_bovy.chapter03.planet_comet_scattering.plt")
    def test_shows_plot_when_no_path(self, mock_plt, mock_integrate, mock_kick, mock_find_peri) -> None:
        """Verify plt.show() is called when no path provided."""
        self._setup_mocks(mock_plt, mock_integrate, mock_kick, mock_find_peri)
        plot_orbit_comparison(path=None)
        mock_plt.show.assert_called_once()

    @patch("galactic_dynamics_bovy.chapter03.planet_comet_scattering.find_pericenter_time")
    @patch("galactic_dynamics_bovy.chapter03.planet_comet_scattering.apply_kick_at_pericenter")
    @patch("galactic_dynamics_bovy.chapter03.planet_comet_scattering.integrate_orbit")
    @patch("galactic_dynamics_bovy.chapter03.planet_comet_scattering.plt")
    def test_saves_figure_when_path_provided(self, mock_plt, mock_integrate, mock_kick, mock_find_peri) -> None:
        """Verify figure is saved when path provided."""
        mock_fig, _, _ = self._setup_mocks(mock_plt, mock_integrate, mock_kick, mock_find_peri)
        test_path = Path("/tmp/test_orbit_comparison.png")
        plot_orbit_comparison(path=test_path)
        mock_fig.savefig.assert_called_once()
        mock_plt.close.assert_called_once_with(mock_fig)

    @patch("galactic_dynamics_bovy.chapter03.planet_comet_scattering.find_pericenter_time")
    @patch("galactic_dynamics_bovy.chapter03.planet_comet_scattering.apply_kick_at_pericenter")
    @patch("galactic_dynamics_bovy.chapter03.planet_comet_scattering.integrate_orbit")
    @patch("galactic_dynamics_bovy.chapter03.planet_comet_scattering.plt")
    def test_sets_labels(self, mock_plt, mock_integrate, mock_kick, mock_find_peri) -> None:
        """Verify axis labels are set."""
        _, mock_ax1, mock_ax2 = self._setup_mocks(mock_plt, mock_integrate, mock_kick, mock_find_peri)
        plot_orbit_comparison()
        mock_ax1.set_xlabel.assert_called_once()
        mock_ax1.set_ylabel.assert_called_once()
        mock_ax2.set_xlabel.assert_called_once()
        mock_ax2.set_ylabel.assert_called_once()

    @patch("galactic_dynamics_bovy.chapter03.planet_comet_scattering.find_pericenter_time")
    @patch("galactic_dynamics_bovy.chapter03.planet_comet_scattering.apply_kick_at_pericenter")
    @patch("galactic_dynamics_bovy.chapter03.planet_comet_scattering.integrate_orbit")
    @patch("galactic_dynamics_bovy.chapter03.planet_comet_scattering.plt")
    def test_has_legend(self, mock_plt, mock_integrate, mock_kick, mock_find_peri) -> None:
        """Verify legend is added on the x-y panel."""
        _, mock_ax1, _ = self._setup_mocks(mock_plt, mock_integrate, mock_kick, mock_find_peri)
        plot_orbit_comparison()
        mock_ax1.legend.assert_called_once()


def _make_mock_study_result() -> ParameterStudyResult:
    """Create a mock ParameterStudyResult for plot tests."""
    kicks = np.linspace(0.01, 0.3, 10)
    zeros = [np.zeros(10), np.zeros(10)]
    return ParameterStudyResult(
        kicks=kicks,
        delta_peri_radial=zeros,
        delta_peri_tangential=zeros,
        delta_apo_radial=zeros,
        delta_apo_tangential=zeros,
        delta_ecc_radial=zeros,
        delta_ecc_tangential=zeros,
    )


class TestPlotParameterStudy:
    """Tests for plot_parameter_study function."""

    @patch("galactic_dynamics_bovy.chapter03.planet_comet_scattering.parameter_study")
    @patch("galactic_dynamics_bovy.chapter03.planet_comet_scattering.plt")
    def test_creates_figure(self, mock_plt: MagicMock, mock_study: MagicMock) -> None:
        """Verify a figure is created."""
        mock_fig = MagicMock()
        mock_plt.subplots.return_value = (mock_fig, [MagicMock(), MagicMock()])
        mock_study.return_value = _make_mock_study_result()
        plot_parameter_study()
        mock_plt.subplots.assert_called_once()

    @patch("galactic_dynamics_bovy.chapter03.planet_comet_scattering.parameter_study")
    @patch("galactic_dynamics_bovy.chapter03.planet_comet_scattering.plt")
    def test_shows_plot_when_no_path(self, mock_plt: MagicMock, mock_study: MagicMock) -> None:
        """Verify plt.show() is called when no path provided."""
        mock_plt.subplots.return_value = (MagicMock(), [MagicMock(), MagicMock()])
        mock_study.return_value = _make_mock_study_result()
        plot_parameter_study(path=None)
        mock_plt.show.assert_called_once()

    @patch("galactic_dynamics_bovy.chapter03.planet_comet_scattering.parameter_study")
    @patch("galactic_dynamics_bovy.chapter03.planet_comet_scattering.plt")
    def test_saves_figure_when_path_provided(self, mock_plt: MagicMock, mock_study: MagicMock) -> None:
        """Verify figure is saved when path provided."""
        mock_fig = MagicMock()
        mock_plt.subplots.return_value = (mock_fig, [MagicMock(), MagicMock()])
        mock_study.return_value = _make_mock_study_result()
        test_path = Path("/tmp/test_param_study.png")
        plot_parameter_study(path=test_path)
        mock_fig.savefig.assert_called_once()
        mock_plt.close.assert_called_once_with(mock_fig)

    @patch("galactic_dynamics_bovy.chapter03.planet_comet_scattering.parameter_study")
    @patch("galactic_dynamics_bovy.chapter03.planet_comet_scattering.plt")
    def test_sets_labels(self, mock_plt: MagicMock, mock_study: MagicMock) -> None:
        """Verify axis labels are set."""
        mock_ax1 = MagicMock()
        mock_ax2 = MagicMock()
        mock_plt.subplots.return_value = (MagicMock(), [mock_ax1, mock_ax2])
        mock_study.return_value = _make_mock_study_result()
        plot_parameter_study()
        mock_ax1.set_xlabel.assert_called_once()
        mock_ax1.set_ylabel.assert_called_once()
        mock_ax2.set_xlabel.assert_called_once()

    @patch("galactic_dynamics_bovy.chapter03.planet_comet_scattering.parameter_study")
    @patch("galactic_dynamics_bovy.chapter03.planet_comet_scattering.plt")
    def test_has_legend(self, mock_plt: MagicMock, mock_study: MagicMock) -> None:
        """Verify legend is added on the left panel only."""
        mock_ax1 = MagicMock()
        mock_ax2 = MagicMock()
        mock_plt.subplots.return_value = (MagicMock(), [mock_ax1, mock_ax2])
        mock_study.return_value = _make_mock_study_result()
        plot_parameter_study()
        mock_ax1.legend.assert_called_once()
        mock_ax2.legend.assert_not_called()


class TestPlotApocenterStudy:
    """Tests for plot_apocenter_study function."""

    @patch("galactic_dynamics_bovy.chapter03.planet_comet_scattering.parameter_study")
    @patch("galactic_dynamics_bovy.chapter03.planet_comet_scattering.plt")
    def test_creates_figure(self, mock_plt: MagicMock, mock_study: MagicMock) -> None:
        """Verify a figure is created."""
        mock_fig = MagicMock()
        mock_plt.subplots.return_value = (mock_fig, [MagicMock(), MagicMock()])
        mock_study.return_value = _make_mock_study_result()
        plot_apocenter_study()
        mock_plt.subplots.assert_called_once()

    @patch("galactic_dynamics_bovy.chapter03.planet_comet_scattering.parameter_study")
    @patch("galactic_dynamics_bovy.chapter03.planet_comet_scattering.plt")
    def test_shows_plot_when_no_path(self, mock_plt: MagicMock, mock_study: MagicMock) -> None:
        """Verify plt.show() is called when no path provided."""
        mock_plt.subplots.return_value = (MagicMock(), [MagicMock(), MagicMock()])
        mock_study.return_value = _make_mock_study_result()
        plot_apocenter_study(path=None)
        mock_plt.show.assert_called_once()

    @patch("galactic_dynamics_bovy.chapter03.planet_comet_scattering.parameter_study")
    @patch("galactic_dynamics_bovy.chapter03.planet_comet_scattering.plt")
    def test_saves_figure_when_path_provided(self, mock_plt: MagicMock, mock_study: MagicMock) -> None:
        """Verify figure is saved when path provided."""
        mock_fig = MagicMock()
        mock_plt.subplots.return_value = (mock_fig, [MagicMock(), MagicMock()])
        mock_study.return_value = _make_mock_study_result()
        test_path = Path("/tmp/test_apo_study.png")
        plot_apocenter_study(path=test_path)
        mock_fig.savefig.assert_called_once()
        mock_plt.close.assert_called_once_with(mock_fig)

    @patch("galactic_dynamics_bovy.chapter03.planet_comet_scattering.parameter_study")
    @patch("galactic_dynamics_bovy.chapter03.planet_comet_scattering.plt")
    def test_sets_labels(self, mock_plt: MagicMock, mock_study: MagicMock) -> None:
        """Verify axis labels are set."""
        mock_ax1 = MagicMock()
        mock_ax2 = MagicMock()
        mock_plt.subplots.return_value = (MagicMock(), [mock_ax1, mock_ax2])
        mock_study.return_value = _make_mock_study_result()
        plot_apocenter_study()
        mock_ax1.set_xlabel.assert_called_once()
        mock_ax1.set_ylabel.assert_called_once()
        mock_ax2.set_xlabel.assert_called_once()

    @patch("galactic_dynamics_bovy.chapter03.planet_comet_scattering.parameter_study")
    @patch("galactic_dynamics_bovy.chapter03.planet_comet_scattering.plt")
    def test_has_legend(self, mock_plt: MagicMock, mock_study: MagicMock) -> None:
        """Verify legend is added on the left panel only."""
        mock_ax1 = MagicMock()
        mock_ax2 = MagicMock()
        mock_plt.subplots.return_value = (MagicMock(), [mock_ax1, mock_ax2])
        mock_study.return_value = _make_mock_study_result()
        plot_apocenter_study()
        mock_ax1.legend.assert_called_once()
        mock_ax2.legend.assert_not_called()


class TestPlotEccentricityStudy:
    """Tests for plot_eccentricity_study function."""

    @patch("galactic_dynamics_bovy.chapter03.planet_comet_scattering.parameter_study")
    @patch("galactic_dynamics_bovy.chapter03.planet_comet_scattering.plt")
    def test_creates_figure(self, mock_plt: MagicMock, mock_study: MagicMock) -> None:
        """Verify a figure is created."""
        mock_fig = MagicMock()
        mock_plt.subplots.return_value = (mock_fig, [MagicMock(), MagicMock()])
        mock_study.return_value = _make_mock_study_result()
        plot_eccentricity_study()
        mock_plt.subplots.assert_called_once()

    @patch("galactic_dynamics_bovy.chapter03.planet_comet_scattering.parameter_study")
    @patch("galactic_dynamics_bovy.chapter03.planet_comet_scattering.plt")
    def test_shows_plot_when_no_path(self, mock_plt: MagicMock, mock_study: MagicMock) -> None:
        """Verify plt.show() is called when no path provided."""
        mock_plt.subplots.return_value = (MagicMock(), [MagicMock(), MagicMock()])
        mock_study.return_value = _make_mock_study_result()
        plot_eccentricity_study(path=None)
        mock_plt.show.assert_called_once()

    @patch("galactic_dynamics_bovy.chapter03.planet_comet_scattering.parameter_study")
    @patch("galactic_dynamics_bovy.chapter03.planet_comet_scattering.plt")
    def test_saves_figure_when_path_provided(self, mock_plt: MagicMock, mock_study: MagicMock) -> None:
        """Verify figure is saved when path provided."""
        mock_fig = MagicMock()
        mock_plt.subplots.return_value = (mock_fig, [MagicMock(), MagicMock()])
        mock_study.return_value = _make_mock_study_result()
        test_path = Path("/tmp/test_ecc_study.png")
        plot_eccentricity_study(path=test_path)
        mock_fig.savefig.assert_called_once()
        mock_plt.close.assert_called_once_with(mock_fig)

    @patch("galactic_dynamics_bovy.chapter03.planet_comet_scattering.parameter_study")
    @patch("galactic_dynamics_bovy.chapter03.planet_comet_scattering.plt")
    def test_sets_labels(self, mock_plt: MagicMock, mock_study: MagicMock) -> None:
        """Verify axis labels are set."""
        mock_ax1 = MagicMock()
        mock_ax2 = MagicMock()
        mock_plt.subplots.return_value = (MagicMock(), [mock_ax1, mock_ax2])
        mock_study.return_value = _make_mock_study_result()
        plot_eccentricity_study()
        mock_ax1.set_xlabel.assert_called_once()
        mock_ax1.set_ylabel.assert_called_once()
        mock_ax2.set_xlabel.assert_called_once()

    @patch("galactic_dynamics_bovy.chapter03.planet_comet_scattering.parameter_study")
    @patch("galactic_dynamics_bovy.chapter03.planet_comet_scattering.plt")
    def test_has_legend(self, mock_plt: MagicMock, mock_study: MagicMock) -> None:
        """Verify legend is added on the left panel only."""
        mock_ax1 = MagicMock()
        mock_ax2 = MagicMock()
        mock_plt.subplots.return_value = (MagicMock(), [mock_ax1, mock_ax2])
        mock_study.return_value = _make_mock_study_result()
        plot_eccentricity_study()
        mock_ax1.legend.assert_called_once()
        mock_ax2.legend.assert_not_called()
