"""Unit tests for IC 2574 rotation curve module."""

from pathlib import Path
from unittest.mock import MagicMock, patch

import numpy as np

from galactic_dynamics_bovy.chapter02.ic_2574 import (
    fit_linear_rotation,
    load_ic2574_data,
    plot_rotation_curve,
)


class TestLoadIc2574Data:
    """Tests for load_ic2574_data function."""

    def test_returns_three_arrays(self) -> None:
        """Verify the function returns radius, velocity, and error arrays."""
        radius, velocity, error = load_ic2574_data()

        assert isinstance(radius, np.ndarray)
        assert isinstance(velocity, np.ndarray)
        assert isinstance(error, np.ndarray)

    def test_arrays_have_same_length(self) -> None:
        """Verify all returned arrays have the same length."""
        radius, velocity, error = load_ic2574_data()

        assert len(radius) == len(velocity) == len(error)

    def test_data_has_expected_size(self) -> None:
        """Verify data has expected number of points."""
        radius, _, _ = load_ic2574_data()

        assert len(radius) == 34

    def test_radius_is_positive(self) -> None:
        """Verify all radius values are positive."""
        radius, _, _ = load_ic2574_data()

        assert np.all(radius > 0)

    def test_velocity_is_positive(self) -> None:
        """Verify all velocity values are positive."""
        _, velocity, _ = load_ic2574_data()

        assert np.all(velocity > 0)

    def test_error_is_non_negative(self) -> None:
        """Verify all error values are non-negative."""
        _, _, error = load_ic2574_data()

        assert np.all(error >= 0)

    def test_radius_range(self) -> None:
        """Verify radius spans expected range."""
        radius, _, _ = load_ic2574_data()

        assert radius.min() < 1.0
        assert radius.max() > 10.0


class TestFitLinearRotation:
    """Tests for fit_linear_rotation function."""

    def test_returns_slope_and_intercept(self) -> None:
        """Verify the function returns slope and intercept."""
        radius = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
        velocity = np.array([10.0, 20.0, 30.0, 40.0, 50.0])

        result = fit_linear_rotation(radius, velocity, r_max=6.0)

        assert len(result) == 2

    def test_perfect_linear_data(self) -> None:
        """Verify fit is exact for perfect linear data."""
        radius = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
        velocity = np.array([15.0, 25.0, 35.0, 45.0, 55.0])  # v = 10*r + 5

        slope, intercept = fit_linear_rotation(radius, velocity, r_max=6.0)

        assert np.isclose(slope, 10.0)
        assert np.isclose(intercept, 5.0)

    def test_r_max_filter(self) -> None:
        """Verify r_max correctly filters data points."""
        radius = np.array([1.0, 2.0, 3.0, 10.0, 20.0])
        velocity = np.array([10.0, 20.0, 30.0, 100.0, 200.0])  # Linear for r < 6

        slope, intercept = fit_linear_rotation(radius, velocity, r_max=6.0)

        # Only points with r < 6 should be used
        assert np.isclose(slope, 10.0)
        assert np.isclose(intercept, 0.0)

    def test_with_real_data(self) -> None:
        """Verify fit works with actual IC 2574 data."""
        radius, velocity, _ = load_ic2574_data()

        slope, intercept = fit_linear_rotation(radius, velocity, r_max=6.0)

        # Expected values from the actual fit
        assert 6.0 < slope < 8.0
        assert 5.0 < intercept < 12.0


class TestPlotRotationCurve:
    """Tests for plot_rotation_curve function."""

    @patch("galactic_dynamics_bovy.chapter02.ic_2574.plt")
    def test_creates_figure(self, mock_plt: MagicMock) -> None:
        """Verify a figure is created."""
        mock_fig = MagicMock()
        mock_ax = MagicMock()
        mock_plt.subplots.return_value = (mock_fig, mock_ax)

        plot_rotation_curve()

        mock_plt.subplots.assert_called_once_with(1, 1, figsize=(6.5, 5))

    @patch("galactic_dynamics_bovy.chapter02.ic_2574.plt")
    def test_shows_plot_when_no_path(self, mock_plt: MagicMock) -> None:
        """Verify plt.show() is called when no path is provided."""
        mock_fig = MagicMock()
        mock_ax = MagicMock()
        mock_plt.subplots.return_value = (mock_fig, mock_ax)

        plot_rotation_curve(path=None)

        mock_plt.show.assert_called_once()

    @patch("galactic_dynamics_bovy.chapter02.ic_2574.plt")
    def test_saves_figure_when_path_provided(self, mock_plt: MagicMock) -> None:
        """Verify figure is saved when path is provided."""
        mock_fig = MagicMock()
        mock_ax = MagicMock()
        mock_plt.subplots.return_value = (mock_fig, mock_ax)
        test_path = Path("/tmp/test_plot.png")

        plot_rotation_curve(path=test_path)

        mock_fig.savefig.assert_called_once()
        mock_plt.close.assert_called_once_with(mock_fig)

    @patch("galactic_dynamics_bovy.chapter02.ic_2574.plt")
    def test_plots_errorbar(self, mock_plt: MagicMock) -> None:
        """Verify errorbar is plotted for observed data."""
        mock_fig = MagicMock()
        mock_ax = MagicMock()
        mock_plt.subplots.return_value = (mock_fig, mock_ax)

        plot_rotation_curve()

        mock_ax.errorbar.assert_called_once()

    @patch("galactic_dynamics_bovy.chapter02.ic_2574.plt")
    def test_plots_linear_fit(self, mock_plt: MagicMock) -> None:
        """Verify linear fit line is plotted."""
        mock_fig = MagicMock()
        mock_ax = MagicMock()
        mock_plt.subplots.return_value = (mock_fig, mock_ax)

        plot_rotation_curve()

        mock_ax.plot.assert_called_once()

    @patch("galactic_dynamics_bovy.chapter02.ic_2574.plt")
    def test_sets_labels(self, mock_plt: MagicMock) -> None:
        """Verify axis labels are set."""
        mock_fig = MagicMock()
        mock_ax = MagicMock()
        mock_plt.subplots.return_value = (mock_fig, mock_ax)

        plot_rotation_curve()

        mock_ax.set_xlabel.assert_called_once()
        mock_ax.set_ylabel.assert_called_once()

    @patch("galactic_dynamics_bovy.chapter02.ic_2574.plt")
    def test_sets_axis_limits(self, mock_plt: MagicMock) -> None:
        """Verify axis limits are set."""
        mock_fig = MagicMock()
        mock_ax = MagicMock()
        mock_plt.subplots.return_value = (mock_fig, mock_ax)

        plot_rotation_curve()

        mock_ax.set_xlim.assert_called_once_with(0, 11)
        mock_ax.set_ylim.assert_called_once_with(0, 80)
