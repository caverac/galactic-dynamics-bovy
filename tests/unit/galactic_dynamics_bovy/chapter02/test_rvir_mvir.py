"""Unit tests for NFW virial mass/radius module."""

from pathlib import Path
from unittest.mock import MagicMock, patch

import numpy as np

from galactic_dynamics_bovy.chapter02.rvir_mvir import (
    _nfw_f,
    get_rvir_mvir_delta,
    plot_rvir_mvir_delta,
)
from galactic_dynamics_bovy.utils.units import GalacticUnits


class TestNfwF:
    """Tests for _nfw_f function."""

    def test_f_at_zero_is_zero(self) -> None:
        """f(0) should be 0."""
        result = _nfw_f(np.array([0.0]))
        assert np.isclose(result[0], 0.0, atol=1e-10)

    def test_f_at_one(self) -> None:
        """f(1) = ln(2) - 0.5 ≈ 0.193."""
        result = _nfw_f(np.array([1.0]))
        expected = np.log(2) - 0.5
        assert np.isclose(result[0], expected)

    def test_f_is_monotonic(self) -> None:
        """f(c) should be monotonically increasing."""
        c = np.linspace(0.1, 100, 100)
        f = _nfw_f(c)
        assert np.all(np.diff(f) > 0)

    def test_f_is_positive_for_positive_c(self) -> None:
        """f(c) > 0 for c > 0."""
        c = np.array([0.1, 1.0, 10.0, 100.0])
        f = _nfw_f(c)
        assert np.all(f > 0)

    def test_f_vectorized(self) -> None:
        """f should work on arrays."""
        c = np.array([1.0, 2.0, 3.0])
        f = _nfw_f(c)
        assert f.shape == c.shape


class TestGetRvirMvirDelta:
    """Tests for get_rvir_mvir_delta function."""

    def test_returns_three_arrays(self) -> None:
        """Function should return (rvir, mvir, delta) arrays."""
        units = GalacticUnits(h=1.0)
        rvir, mvir, delta = get_rvir_mvir_delta(0.00035, 16.0, units.rho_crit)

        assert isinstance(rvir, np.ndarray)
        assert isinstance(mvir, np.ndarray)
        assert isinstance(delta, np.ndarray)

    def test_arrays_have_same_length(self) -> None:
        """All returned arrays should have same length."""
        units = GalacticUnits(h=1.0)
        rvir, mvir, delta = get_rvir_mvir_delta(0.00035, 16.0, units.rho_crit, n_points=100)

        assert len(rvir) == len(mvir) == len(delta) == 100

    def test_rvir_range(self) -> None:
        """rvir should span the requested range."""
        units = GalacticUnits(h=1.0)
        rvir, _, _ = get_rvir_mvir_delta(0.00035, 16.0, units.rho_crit, rvir_range=(100, 200))

        assert np.isclose(rvir.min(), 100)
        assert np.isclose(rvir.max(), 200)

    def test_mvir_is_positive(self) -> None:
        """Virial mass should be positive."""
        units = GalacticUnits(h=1.0)
        _, mvir, _ = get_rvir_mvir_delta(0.00035, 16.0, units.rho_crit)

        assert np.all(mvir > 0)

    def test_delta_is_positive(self) -> None:
        """Overdensity should be positive."""
        units = GalacticUnits(h=1.0)
        _, _, delta = get_rvir_mvir_delta(0.00035, 16.0, units.rho_crit)

        assert np.all(delta > 0)

    def test_delta_decreases_with_rvir(self) -> None:
        """Larger rvir should give smaller overdensity."""
        units = GalacticUnits(h=1.0)
        _, _, delta = get_rvir_mvir_delta(0.00035, 16.0, units.rho_crit)

        # Since rvir is increasing, delta should be decreasing
        assert np.all(np.diff(delta) < 0)

    def test_mvir_increases_with_rvir(self) -> None:
        """Larger rvir should enclose more mass."""
        units = GalacticUnits(h=1.0)
        _, mvir, _ = get_rvir_mvir_delta(0.00035, 16.0, units.rho_crit)

        assert np.all(np.diff(mvir) > 0)

    def test_delta_200_gives_reasonable_rvir(self) -> None:
        """At Δ≈200, rvir should be around 90-100 kpc for MW-like halo."""
        units = GalacticUnits(h=1.0)
        rvir, _, delta = get_rvir_mvir_delta(0.00035, 16.0, units.rho_crit, rvir_range=(50, 200))

        # Find where delta is closest to 200
        idx = np.argmin(np.abs(delta - 200))
        rvir_200 = rvir[idx]

        assert 80 < rvir_200 < 120

    def test_mvir_formula_consistency(self) -> None:
        """M_vir = (4π/3) r³ ρ_crit Δ should hold."""
        units = GalacticUnits(h=1.0)
        rvir, mvir, delta = get_rvir_mvir_delta(0.00035, 16.0, units.rho_crit)

        mvir_check = 4 * np.pi * rvir**3 * units.rho_crit * delta / 3
        assert np.allclose(mvir, mvir_check)


class TestPlotRvirMvirDelta:
    """Tests for plot_rvir_mvir_delta function."""

    @patch("galactic_dynamics_bovy.chapter02.rvir_mvir.plt")
    def test_creates_figure(self, mock_plt: MagicMock) -> None:
        """Verify a figure is created with correct dimensions."""
        mock_fig = MagicMock()
        mock_ax1 = MagicMock()
        mock_ax2 = MagicMock()
        mock_plt.subplots.return_value = (mock_fig, [mock_ax1, mock_ax2])

        plot_rvir_mvir_delta()

        mock_plt.subplots.assert_called_once_with(2, 1, figsize=(6.5, 7), sharex=True)

    @patch("galactic_dynamics_bovy.chapter02.rvir_mvir.plt")
    def test_shows_plot_when_no_path(self, mock_plt: MagicMock) -> None:
        """Verify plt.show() is called when no path provided."""
        mock_fig = MagicMock()
        mock_ax1 = MagicMock()
        mock_ax2 = MagicMock()
        mock_plt.subplots.return_value = (mock_fig, [mock_ax1, mock_ax2])

        plot_rvir_mvir_delta(path=None)

        mock_plt.show.assert_called_once()

    @patch("galactic_dynamics_bovy.chapter02.rvir_mvir.plt")
    def test_saves_figure_when_path_provided(self, mock_plt: MagicMock) -> None:
        """Verify figure is saved when path provided."""
        mock_fig = MagicMock()
        mock_ax1 = MagicMock()
        mock_ax2 = MagicMock()
        mock_plt.subplots.return_value = (mock_fig, [mock_ax1, mock_ax2])
        test_path = Path("/tmp/test_plot.png")

        plot_rvir_mvir_delta(path=test_path)

        mock_fig.savefig.assert_called_once()
        mock_plt.close.assert_called_once_with(mock_fig)

    @patch("galactic_dynamics_bovy.chapter02.rvir_mvir.plt")
    def test_plots_two_curves(self, mock_plt: MagicMock) -> None:
        """Verify both panels have a plot call."""
        mock_fig = MagicMock()
        mock_ax1 = MagicMock()
        mock_ax2 = MagicMock()
        mock_plt.subplots.return_value = (mock_fig, [mock_ax1, mock_ax2])

        plot_rvir_mvir_delta()

        mock_ax1.plot.assert_called_once()
        mock_ax2.plot.assert_called_once()

    @patch("galactic_dynamics_bovy.chapter02.rvir_mvir.plt")
    def test_sets_labels(self, mock_plt: MagicMock) -> None:
        """Verify axis labels are set."""
        mock_fig = MagicMock()
        mock_ax1 = MagicMock()
        mock_ax2 = MagicMock()
        mock_plt.subplots.return_value = (mock_fig, [mock_ax1, mock_ax2])

        plot_rvir_mvir_delta()

        # Top panel: y-label only (shared x-axis)
        mock_ax1.set_ylabel.assert_called_once()
        # Bottom panel: both labels
        mock_ax2.set_xlabel.assert_called_once()
        mock_ax2.set_ylabel.assert_called_once()

    @patch("galactic_dynamics_bovy.chapter02.rvir_mvir.plt")
    def test_enables_grid(self, mock_plt: MagicMock) -> None:
        """Verify grid is enabled on both panels."""
        mock_fig = MagicMock()
        mock_ax1 = MagicMock()
        mock_ax2 = MagicMock()
        mock_plt.subplots.return_value = (mock_fig, [mock_ax1, mock_ax2])

        plot_rvir_mvir_delta()

        mock_ax1.grid.assert_called_once()
        mock_ax2.grid.assert_called_once()
