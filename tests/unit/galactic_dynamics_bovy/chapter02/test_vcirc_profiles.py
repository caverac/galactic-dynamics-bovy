"""Unit tests for circular velocity profiles module."""

from pathlib import Path
from unittest.mock import MagicMock, patch

import numpy as np

from galactic_dynamics_bovy.chapter02.vcirc_profiles import (
    plot_vcirc_profiles,
    vcirc_hernquist,
    vcirc_jaffe,
    vcirc_nfw,
)


class TestVcircJaffe:
    """Tests for vcirc_jaffe function."""

    def test_at_zero(self) -> None:
        """v_c(0) = 1 for Jaffe."""
        x = np.array([0.0])
        result = vcirc_jaffe(x)
        assert np.isclose(result[0], 1.0)

    def test_at_one(self) -> None:
        """v_c(1) = 1/sqrt(2) for Jaffe."""
        x = np.array([1.0])
        result = vcirc_jaffe(x)
        assert np.isclose(result[0], 1.0 / np.sqrt(2.0))

    def test_decreasing(self) -> None:
        """Jaffe rotation curve is monotonically decreasing."""
        x = np.logspace(-2, 2, 100)
        v = vcirc_jaffe(x)
        assert np.all(np.diff(v) < 0)

    def test_vectorized(self) -> None:
        """Function should work on arrays."""
        x = np.array([0.1, 1.0, 10.0])
        v = vcirc_jaffe(x)
        assert v.shape == x.shape


class TestVcircHernquist:
    """Tests for vcirc_hernquist function."""

    def test_at_zero(self) -> None:
        """v_c(0) = 0 for Hernquist."""
        x = np.array([1e-10])
        result = vcirc_hernquist(x)
        assert result[0] < 0.01

    def test_at_one(self) -> None:
        """v_c(1) = 1/(2*sqrt(2)) for Hernquist."""
        x = np.array([1.0])
        result = vcirc_hernquist(x)
        expected = np.sqrt(1.0 / (2.0 * 4.0))  # x/(2*(1+x)^2) at x=1
        assert np.isclose(result[0], expected)

    def test_has_peak(self) -> None:
        """Hernquist rotation curve has a peak."""
        x = np.logspace(-2, 2, 1000)
        v = vcirc_hernquist(x)
        # Find max is not at edges
        idx_max = np.argmax(v)
        assert 0 < idx_max < len(x) - 1

    def test_vectorized(self) -> None:
        """Function should work on arrays."""
        x = np.array([0.1, 1.0, 10.0])
        v = vcirc_hernquist(x)
        assert v.shape == x.shape


class TestVcircNfw:
    """Tests for vcirc_nfw function."""

    def test_at_small_x(self) -> None:
        """v_c approaches 0 as x -> 0 for NFW."""
        x = np.array([1e-5])
        result = vcirc_nfw(x)
        assert result[0] < 0.01

    def test_at_one(self) -> None:
        """v_c(1) uses f(1) = ln(2) - 0.5."""
        x = np.array([1.0])
        result = vcirc_nfw(x)
        f_1 = np.log(2) - 0.5
        expected = np.sqrt(f_1)
        assert np.isclose(result[0], expected)

    def test_has_peak(self) -> None:
        """NFW rotation curve has a peak."""
        x = np.logspace(-2, 2, 1000)
        v = vcirc_nfw(x)
        # Find max is not at edges
        idx_max = np.argmax(v)
        assert 0 < idx_max < len(x) - 1

    def test_vectorized(self) -> None:
        """Function should work on arrays."""
        x = np.array([0.1, 1.0, 10.0])
        v = vcirc_nfw(x)
        assert v.shape == x.shape


class TestPlotVcircProfiles:
    """Tests for plot_vcirc_profiles function."""

    @patch("galactic_dynamics_bovy.chapter02.vcirc_profiles.plt")
    def test_creates_figure(self, mock_plt: MagicMock) -> None:
        """Verify a figure is created."""
        mock_fig = MagicMock()
        mock_ax = MagicMock()
        mock_plt.subplots.return_value = (mock_fig, mock_ax)

        plot_vcirc_profiles()

        mock_plt.subplots.assert_called_once()

    @patch("galactic_dynamics_bovy.chapter02.vcirc_profiles.plt")
    def test_shows_plot_when_no_path(self, mock_plt: MagicMock) -> None:
        """Verify plt.show() is called when no path provided."""
        mock_fig = MagicMock()
        mock_ax = MagicMock()
        mock_plt.subplots.return_value = (mock_fig, mock_ax)

        plot_vcirc_profiles(path=None)

        mock_plt.show.assert_called_once()

    @patch("galactic_dynamics_bovy.chapter02.vcirc_profiles.plt")
    def test_saves_figure_when_path_provided(self, mock_plt: MagicMock) -> None:
        """Verify figure is saved when path provided."""
        mock_fig = MagicMock()
        mock_ax = MagicMock()
        mock_plt.subplots.return_value = (mock_fig, mock_ax)
        test_path = Path("/tmp/test_plot.png")

        plot_vcirc_profiles(path=test_path)

        mock_fig.savefig.assert_called_once()
        mock_plt.close.assert_called_once_with(mock_fig)

    @patch("galactic_dynamics_bovy.chapter02.vcirc_profiles.plt")
    def test_plots_three_curves(self, mock_plt: MagicMock) -> None:
        """Verify all three profiles are plotted."""
        mock_fig = MagicMock()
        mock_ax = MagicMock()
        mock_plt.subplots.return_value = (mock_fig, mock_ax)

        plot_vcirc_profiles()

        assert mock_ax.plot.call_count == 3

    @patch("galactic_dynamics_bovy.chapter02.vcirc_profiles.plt")
    def test_sets_labels(self, mock_plt: MagicMock) -> None:
        """Verify axis labels are set."""
        mock_fig = MagicMock()
        mock_ax = MagicMock()
        mock_plt.subplots.return_value = (mock_fig, mock_ax)

        plot_vcirc_profiles()

        mock_ax.set_xlabel.assert_called_once()
        mock_ax.set_ylabel.assert_called_once()

    @patch("galactic_dynamics_bovy.chapter02.vcirc_profiles.plt")
    def test_has_legend(self, mock_plt: MagicMock) -> None:
        """Verify legend is added."""
        mock_fig = MagicMock()
        mock_ax = MagicMock()
        mock_plt.subplots.return_value = (mock_fig, mock_ax)

        plot_vcirc_profiles()

        mock_ax.legend.assert_called_once()
