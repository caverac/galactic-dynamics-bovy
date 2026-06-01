"""Unit tests for the isothermal Gaussian-limit demonstration (Problem 5.9)."""

from pathlib import Path
from unittest.mock import MagicMock, patch

import numpy as np

from galactic_dynamics_bovy.chapter05.isothermal_gaussian_limit import (
    _setup_axes,
    central_excess_kurtosis,
    plot_king_maxwellian_convergence,
    W0_VALUES,
)

MODULE = "galactic_dynamics_bovy.chapter05.isothermal_gaussian_limit"


class TestCentralExcessKurtosis:
    """Tests for central_excess_kurtosis (real galpy: fast, deterministic)."""

    def test_returns_float(self) -> None:
        """Should return a Python float."""
        result = central_excess_kurtosis(4.0)
        assert isinstance(result, float)

    def test_converges_toward_gaussian(self) -> None:
        """Higher W0 yields a central distribution closer to Gaussian."""
        low = central_excess_kurtosis(2.0)
        high = central_excess_kurtosis(20.0)
        # The barely-bound model is strongly platykurtic (negative excess);
        # the near-isothermal model is essentially Gaussian.
        assert low < -0.4
        assert abs(high) < 0.01

    def test_monotonic_approach(self) -> None:
        """Excess kurtosis increases monotonically toward zero with W0."""
        kappa = [central_excess_kurtosis(w0) for w0 in (1.0, 4.0, 8.0, 16.0)]
        assert all(np.diff(kappa) > 0)
        assert all(k < 0.001 for k in kappa)

    def test_maxwellian_identity(self) -> None:
        """A pure exponential DF (no truncation effect) gives near-zero kurtosis."""
        # With a deep well the truncation is negligible and kappa -> 0.
        assert abs(central_excess_kurtosis(20.0)) < 1e-3


class TestSetupAxes:
    """Tests for the _setup_axes helper."""

    def test_applies_formatting(self) -> None:
        """Should configure ticks and grid."""
        mock_ax = MagicMock()
        _setup_axes(mock_ax)
        assert mock_ax.tick_params.called
        mock_ax.grid.assert_called_once()


class TestPlotKingMaxwellianConvergence:
    """Tests for the top-level asset plot function."""

    @patch(f"{MODULE}.central_excess_kurtosis")
    @patch(f"{MODULE}.plt")
    def test_creates_figure(self, mock_plt: MagicMock, mock_kurtosis: MagicMock) -> None:
        """A single-panel figure is created with one kurtosis per W0."""
        mock_plt.subplots.return_value = (MagicMock(), MagicMock())
        mock_kurtosis.side_effect = [-1.0 / w0 for w0 in W0_VALUES]
        plot_king_maxwellian_convergence()
        mock_plt.subplots.assert_called_once()
        assert mock_kurtosis.call_count == len(W0_VALUES)

    @patch(f"{MODULE}.central_excess_kurtosis")
    @patch(f"{MODULE}.plt")
    def test_plots_curve_and_reference(self, mock_plt: MagicMock, mock_kurtosis: MagicMock) -> None:
        """The kurtosis curve and the zero reference line are drawn."""
        mock_ax = MagicMock()
        mock_plt.subplots.return_value = (MagicMock(), mock_ax)
        mock_kurtosis.side_effect = [-1.0 / w0 for w0 in W0_VALUES]
        plot_king_maxwellian_convergence()
        mock_ax.plot.assert_called_once()
        mock_ax.axhline.assert_called_once()
        mock_ax.set_xlabel.assert_called_once()

    @patch(f"{MODULE}.central_excess_kurtosis")
    @patch(f"{MODULE}.plt")
    def test_shows_when_no_path(self, mock_plt: MagicMock, mock_kurtosis: MagicMock) -> None:
        """plt.show() is called when no path is given."""
        mock_plt.subplots.return_value = (MagicMock(), MagicMock())
        mock_kurtosis.side_effect = [-1.0 / w0 for w0 in W0_VALUES]
        plot_king_maxwellian_convergence(path=None)
        mock_plt.show.assert_called_once()

    @patch(f"{MODULE}.save_figure_if_changed")
    @patch(f"{MODULE}.central_excess_kurtosis")
    @patch(f"{MODULE}.plt")
    def test_saves_when_path_provided(
        self, mock_plt: MagicMock, mock_kurtosis: MagicMock, mock_save: MagicMock, tmp_path: Path
    ) -> None:
        """The figure is saved and closed when a path is provided."""
        mock_fig = MagicMock()
        mock_plt.subplots.return_value = (mock_fig, MagicMock())
        mock_kurtosis.side_effect = [-1.0 / w0 for w0 in W0_VALUES]
        plot_king_maxwellian_convergence(path=tmp_path / "fig.png")
        mock_save.assert_called_once()
        mock_plt.close.assert_called_once_with(mock_fig)
