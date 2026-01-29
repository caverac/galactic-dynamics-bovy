"""Unit tests for spherical exponential disk module."""

from pathlib import Path
from unittest.mock import MagicMock, patch

import numpy as np

from galactic_dynamics_bovy.chapter02.spherical_exponential import (
    plot_exponential_vcirc,
    vcirc_exponential_disk,
    vcirc_spherical_exponential,
)


class TestVcircSphericalExponential:
    """Tests for vcirc_spherical_exponential function."""

    def test_at_small_x(self) -> None:
        """v_c^2 approaches 0 as x -> 0."""
        x = np.array([0.01])
        result = vcirc_spherical_exponential(x)
        assert result[0] < 0.01

    def test_at_large_x(self) -> None:
        """v_c^2 approaches 1/x as x -> infinity."""
        x = np.array([100.0])
        result = vcirc_spherical_exponential(x)
        # At large x, (1 - (1+x)e^{-x})/x -> 1/x
        assert np.isclose(result[0], 1.0 / x[0], rtol=0.01)

    def test_has_peak(self) -> None:
        """Rotation curve has a peak."""
        x = np.linspace(0.1, 10, 1000)
        v2 = vcirc_spherical_exponential(x)
        # Find max is not at edges
        idx_max = np.argmax(v2)
        assert 0 < idx_max < len(x) - 1

    def test_peak_location(self) -> None:
        """Peak should be around x ~ 2."""
        x = np.linspace(0.1, 10, 1000)
        v2 = vcirc_spherical_exponential(x)
        x_peak = x[np.argmax(v2)]
        assert 1.5 < x_peak < 2.5

    def test_vectorized(self) -> None:
        """Function should work on arrays."""
        x = np.array([0.5, 1.0, 2.0, 5.0])
        v2 = vcirc_spherical_exponential(x)
        assert v2.shape == x.shape


class TestVcircExponentialDisk:
    """Tests for vcirc_exponential_disk function."""

    def test_at_small_x(self) -> None:
        """v_c^2 approaches 0 as x -> 0."""
        x = np.array([0.01])
        result = vcirc_exponential_disk(x)
        assert result[0] < 0.01

    def test_has_peak(self) -> None:
        """Rotation curve has a peak."""
        x = np.linspace(0.1, 10, 1000)
        v2 = vcirc_exponential_disk(x)
        # Find max is not at edges
        idx_max = np.argmax(v2)
        assert 0 < idx_max < len(x) - 1

    def test_peak_location(self) -> None:
        """Peak should be around x ~ 2.2 for exponential disk."""
        x = np.linspace(0.1, 10, 1000)
        v2 = vcirc_exponential_disk(x)
        x_peak = x[np.argmax(v2)]
        assert 2.0 < x_peak < 2.5

    def test_vectorized(self) -> None:
        """Function should work on arrays."""
        x = np.array([0.5, 1.0, 2.0, 5.0])
        v2 = vcirc_exponential_disk(x)
        assert v2.shape == x.shape

    def test_positive(self) -> None:
        """v_c^2 should be positive for all x > 0."""
        x = np.linspace(0.1, 10, 100)
        v2 = vcirc_exponential_disk(x)
        assert np.all(v2 > 0)


class TestPlotExponentialVcirc:
    """Tests for plot_exponential_vcirc function."""

    @patch("galactic_dynamics_bovy.chapter02.spherical_exponential.plt")
    def test_creates_figure(self, mock_plt: MagicMock) -> None:
        """Verify a figure is created."""
        mock_fig = MagicMock()
        mock_ax = MagicMock()
        mock_plt.subplots.return_value = (mock_fig, mock_ax)

        plot_exponential_vcirc()

        mock_plt.subplots.assert_called_once()

    @patch("galactic_dynamics_bovy.chapter02.spherical_exponential.plt")
    def test_shows_plot_when_no_path(self, mock_plt: MagicMock) -> None:
        """Verify plt.show() is called when no path provided."""
        mock_fig = MagicMock()
        mock_ax = MagicMock()
        mock_plt.subplots.return_value = (mock_fig, mock_ax)

        plot_exponential_vcirc(path=None)

        mock_plt.show.assert_called_once()

    @patch("galactic_dynamics_bovy.chapter02.spherical_exponential.plt")
    def test_saves_figure_when_path_provided(self, mock_plt: MagicMock) -> None:
        """Verify figure is saved when path provided."""
        mock_fig = MagicMock()
        mock_ax = MagicMock()
        mock_plt.subplots.return_value = (mock_fig, mock_ax)
        test_path = Path("/tmp/test_plot.png")

        plot_exponential_vcirc(path=test_path)

        mock_fig.savefig.assert_called_once()
        mock_plt.close.assert_called_once_with(mock_fig)

    @patch("galactic_dynamics_bovy.chapter02.spherical_exponential.plt")
    def test_plots_two_curves(self, mock_plt: MagicMock) -> None:
        """Verify both profiles are plotted."""
        mock_fig = MagicMock()
        mock_ax = MagicMock()
        mock_plt.subplots.return_value = (mock_fig, mock_ax)

        plot_exponential_vcirc()

        assert mock_ax.plot.call_count == 2

    @patch("galactic_dynamics_bovy.chapter02.spherical_exponential.plt")
    def test_sets_labels(self, mock_plt: MagicMock) -> None:
        """Verify axis labels are set."""
        mock_fig = MagicMock()
        mock_ax = MagicMock()
        mock_plt.subplots.return_value = (mock_fig, mock_ax)

        plot_exponential_vcirc()

        mock_ax.set_xlabel.assert_called_once()
        mock_ax.set_ylabel.assert_called_once()

    @patch("galactic_dynamics_bovy.chapter02.spherical_exponential.plt")
    def test_has_legend(self, mock_plt: MagicMock) -> None:
        """Verify legend is added."""
        mock_fig = MagicMock()
        mock_ax = MagicMock()
        mock_plt.subplots.return_value = (mock_fig, mock_ax)

        plot_exponential_vcirc()

        mock_ax.legend.assert_called_once()
