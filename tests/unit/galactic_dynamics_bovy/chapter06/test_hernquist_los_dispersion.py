"""Unit tests for the Hernquist line-of-sight dispersion (Problem 6.7)."""

from pathlib import Path
from unittest.mock import MagicMock, patch

import numpy as np

from galactic_dynamics_bovy.chapter06.hernquist_los_dispersion import (
    _setup_axes,
    BETAS,
    hernquist_mass,
    plot_hernquist_los_dispersion,
    sigma_los_profile,
    surface_density,
)

MODULE = "galactic_dynamics_bovy.chapter06.hernquist_los_dispersion"


class TestHernquistMass:
    """Tests for the Hernquist enclosed mass."""

    def test_limits(self) -> None:
        """M(<r) -> 0 at the center and -> M (=1) at infinity."""
        assert hernquist_mass(1e-6) < 1e-10
        assert np.isclose(hernquist_mass(1e6), 1.0, atol=1e-4)

    def test_half_mass_radius(self) -> None:
        """M(<r) = r^2/(1+r)^2 gives M = 1/4 at r = a."""
        assert np.isclose(hernquist_mass(1.0), 0.25)


class TestSurfaceDensity:
    """Tests for the projected surface density."""

    def test_positive_and_decreasing(self) -> None:
        """Surface density is positive and falls with projected radius."""
        assert surface_density(0.1) > surface_density(1.0) > surface_density(5.0) > 0.0


class TestSigmaLosProfile:
    """Tests for the projected dispersion profile and the anisotropy signature."""

    def test_isotropic_central_value(self) -> None:
        """Isotropic central sigma_los ~ 0.3 sqrt(GM/a) (Hernquist 1990)."""
        sigma = sigma_los_profile(np.array([0.05]), beta=0.0, n_grid=160)
        assert np.isclose(sigma[0], 0.30, atol=0.03)

    def test_radial_hotter_than_tangential_in_center(self) -> None:
        """Radial anisotropy raises, tangential lowers, the central sigma_los."""
        radius = np.array([0.05])
        radial = sigma_los_profile(radius, beta=0.5, n_grid=120)[0]
        iso = sigma_los_profile(radius, beta=0.0, n_grid=120)[0]
        tang = sigma_los_profile(radius, beta=-0.5, n_grid=120)[0]
        assert radial > iso > tang

    def test_positive_and_finite(self) -> None:
        """The profile is positive and finite across radius."""
        sigma = sigma_los_profile(np.array([0.1, 1.0, 5.0]), beta=0.0, n_grid=80)
        assert np.all(sigma > 0.0) and np.all(np.isfinite(sigma))


class TestSetupAxes:
    """Tests for the _setup_axes helper."""

    def test_applies_formatting(self) -> None:
        """Should configure ticks."""
        mock_ax = MagicMock()
        _setup_axes(mock_ax)
        assert mock_ax.tick_params.called


class TestPlotHernquistLosDispersion:
    """Tests for the asset plot function."""

    @patch(f"{MODULE}.sigma_los_profile")
    @patch(f"{MODULE}.plt")
    def test_creates_one_line_per_beta(self, mock_plt: MagicMock, mock_profile: MagicMock) -> None:
        """One curve is drawn per anisotropy value."""
        mock_ax = MagicMock()
        mock_plt.subplots.return_value = (MagicMock(), mock_ax)
        mock_profile.return_value = np.full(40, 0.3)
        plot_hernquist_los_dispersion()
        mock_plt.subplots.assert_called_once()
        assert mock_ax.plot.call_count == len(BETAS)
        assert mock_profile.call_count == len(BETAS)

    @patch(f"{MODULE}.sigma_los_profile")
    @patch(f"{MODULE}.plt")
    def test_shows_when_no_path(self, mock_plt: MagicMock, mock_profile: MagicMock) -> None:
        """plt.show() is called when no path is given."""
        mock_plt.subplots.return_value = (MagicMock(), MagicMock())
        mock_profile.return_value = np.full(40, 0.3)
        plot_hernquist_los_dispersion(path=None)
        mock_plt.show.assert_called_once()

    @patch(f"{MODULE}.save_figure_if_changed")
    @patch(f"{MODULE}.sigma_los_profile")
    @patch(f"{MODULE}.plt")
    def test_saves_when_path_provided(
        self, mock_plt: MagicMock, mock_profile: MagicMock, mock_save: MagicMock, tmp_path: Path
    ) -> None:
        """The figure is saved and closed when a path is provided."""
        mock_fig = MagicMock()
        mock_plt.subplots.return_value = (mock_fig, MagicMock())
        mock_profile.return_value = np.full(40, 0.3)
        plot_hernquist_los_dispersion(path=tmp_path / "fig.png")
        mock_save.assert_called_once()
        mock_plt.close.assert_called_once_with(mock_fig)
