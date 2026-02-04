"""Unit tests for effective radius gamma module."""

from pathlib import Path
from unittest.mock import MagicMock, patch
import warnings

import numpy as np
from scipy.integrate import IntegrationWarning

from galactic_dynamics_bovy.chapter02.effective_radius_gamma import (
    compute_effective_radius,
    compute_s_direct,
    MInterpolator,
    plot_effective_radius_gamma,
    SInterpolator,
)


class TestComputeSDirect:
    """Tests for compute_s_direct function."""

    def test_positive_for_positive_x(self) -> None:
        """s(x) should be positive for positive x."""
        result = compute_s_direct(1.0, gamma=1.0)
        assert result > 0

    def test_decreases_with_x(self) -> None:
        """s(x) should decrease as x increases."""
        s1 = compute_s_direct(0.1, gamma=1.0)
        s2 = compute_s_direct(1.0, gamma=1.0)
        s3 = compute_s_direct(10.0, gamma=1.0)
        assert s1 > s2 > s3

    def test_scales_with_gamma(self) -> None:
        """Different gamma values should give different results."""
        s1 = compute_s_direct(1.0, gamma=0.5)
        s2 = compute_s_direct(1.0, gamma=1.5)
        assert s1 != s2

    def test_handles_x_zero(self) -> None:
        """Should handle x=0 case (may produce integration warning due to singularity)."""
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=IntegrationWarning)
            result = compute_s_direct(0.0, gamma=1.0)
        assert np.isfinite(result)
        assert result > 0


class TestSInterpolator:
    """Tests for SInterpolator class."""

    def test_initialization(self) -> None:
        """Should initialize without error."""
        interp = SInterpolator(gamma=1.0, n_points=50)
        assert interp.gamma == 1.0

    def test_call_scalar(self) -> None:
        """Should return scalar for scalar input."""
        interp = SInterpolator(gamma=1.0, n_points=50)
        result = interp(1.0)
        assert isinstance(result, float)
        assert result > 0

    def test_call_array(self) -> None:
        """Should return array for array input."""
        interp = SInterpolator(gamma=1.0, n_points=50)
        x = np.array([0.1, 1.0, 10.0])
        result = interp(x)
        assert isinstance(result, np.ndarray)
        assert len(result) == 3
        assert all(result > 0)

    def test_matches_direct_in_range(self) -> None:
        """Interpolated values should match direct computation."""
        interp = SInterpolator(gamma=1.0, n_points=100)
        x = 1.0
        interpolated = interp(x)
        direct = compute_s_direct(x, gamma=1.0)
        assert np.isclose(interpolated, direct, rtol=0.01)

    def test_outside_range_uses_direct(self) -> None:
        """Values outside range should use direct computation."""
        interp = SInterpolator(gamma=1.0, x_min=0.1, x_max=10.0, n_points=50)
        # Test value outside range
        result = interp(100.0)
        direct = compute_s_direct(100.0, gamma=1.0)
        assert np.isclose(result, direct, rtol=0.01)


class TestMInterpolator:
    """Tests for MInterpolator class."""

    def test_initialization(self) -> None:
        """Should initialize without error."""
        s_interp = SInterpolator(gamma=1.0, n_points=50)
        m_interp = MInterpolator(s_interp)
        assert m_interp.s_interp is s_interp

    def test_call_scalar(self) -> None:
        """Should return scalar for scalar input."""
        s_interp = SInterpolator(gamma=1.0, n_points=50)
        m_interp = MInterpolator(s_interp)
        result = m_interp(1.0)
        assert isinstance(result, float)

    def test_monotonically_increasing(self) -> None:
        """m(x) should be monotonically increasing."""
        s_interp = SInterpolator(gamma=1.0, n_points=50)
        m_interp = MInterpolator(s_interp)
        x = np.array([0.01, 0.1, 1.0, 10.0])
        m = m_interp(x)
        assert all(np.diff(m) > 0)

    def test_inverse(self) -> None:
        """Inverse should return x such that m(x) = target."""
        s_interp = SInterpolator(gamma=1.0, n_points=100)
        m_interp = MInterpolator(s_interp)
        # Get m at a known x
        x_test = 1.0
        m_test = float(m_interp(x_test))
        # Invert
        x_recovered = m_interp.inverse(m_test)
        assert np.isclose(x_recovered, x_test, rtol=0.05)


class TestComputeEffectiveRadius:
    """Tests for compute_effective_radius function."""

    def test_returns_positive(self) -> None:
        """Effective radius should be positive."""
        xe = compute_effective_radius(gamma=1.0)
        assert xe > 0

    def test_returns_finite(self) -> None:
        """Effective radius should be finite."""
        xe = compute_effective_radius(gamma=1.0)
        assert np.isfinite(xe)

    def test_varies_with_gamma(self) -> None:
        """Different gamma values should give different effective radii."""
        xe1 = compute_effective_radius(gamma=0.5)
        xe2 = compute_effective_radius(gamma=1.5)
        assert xe1 != xe2

    def test_unnormalized_returns_positive(self) -> None:
        """Unnormalized effective radius should be positive."""
        xe = compute_effective_radius(gamma=1.0, normalize=False)
        assert xe > 0
        assert np.isfinite(xe)


class TestPlotEffectiveRadiusGamma:
    """Tests for plot_effective_radius_gamma function."""

    @patch("galactic_dynamics_bovy.chapter02.effective_radius_gamma.plt")
    @patch("galactic_dynamics_bovy.chapter02.effective_radius_gamma.compute_effective_radius")
    def test_creates_figure(
        self,
        mock_compute: MagicMock,
        mock_plt: MagicMock,
    ) -> None:
        """Verify a figure is created."""
        mock_fig = MagicMock()
        mock_ax = MagicMock()
        mock_plt.subplots.return_value = (mock_fig, mock_ax)
        mock_compute.return_value = 0.5

        plot_effective_radius_gamma()

        mock_plt.subplots.assert_called_once()

    @patch("galactic_dynamics_bovy.chapter02.effective_radius_gamma.plt")
    @patch("galactic_dynamics_bovy.chapter02.effective_radius_gamma.compute_effective_radius")
    def test_shows_plot_when_no_path(
        self,
        mock_compute: MagicMock,
        mock_plt: MagicMock,
    ) -> None:
        """Verify plt.show() is called when no path provided."""
        mock_fig = MagicMock()
        mock_ax = MagicMock()
        mock_plt.subplots.return_value = (mock_fig, mock_ax)
        mock_compute.return_value = 0.5

        plot_effective_radius_gamma(path=None)

        mock_plt.show.assert_called_once()

    @patch("galactic_dynamics_bovy.chapter02.effective_radius_gamma.plt")
    @patch("galactic_dynamics_bovy.chapter02.effective_radius_gamma.save_figure_if_changed")
    @patch("galactic_dynamics_bovy.chapter02.effective_radius_gamma.compute_effective_radius")
    def test_saves_figure_when_path_provided(
        self,
        mock_compute: MagicMock,
        mock_save: MagicMock,
        mock_plt: MagicMock,
    ) -> None:
        """Verify figure is saved when path provided."""
        mock_fig = MagicMock()
        mock_ax = MagicMock()
        mock_plt.subplots.return_value = (mock_fig, mock_ax)
        mock_compute.return_value = 0.5
        test_path = Path("/tmp/test_plot.png")

        plot_effective_radius_gamma(path=test_path)

        mock_save.assert_called_once()
        mock_plt.close.assert_called_once_with(mock_fig)

    @patch("galactic_dynamics_bovy.chapter02.effective_radius_gamma.plt")
    @patch("galactic_dynamics_bovy.chapter02.effective_radius_gamma.compute_effective_radius")
    def test_sets_labels(
        self,
        mock_compute: MagicMock,
        mock_plt: MagicMock,
    ) -> None:
        """Verify axis labels are set."""
        mock_fig = MagicMock()
        mock_ax = MagicMock()
        mock_plt.subplots.return_value = (mock_fig, mock_ax)
        mock_compute.return_value = 0.5

        plot_effective_radius_gamma()

        mock_ax.set_xlabel.assert_called_once()
        mock_ax.set_ylabel.assert_called_once()
