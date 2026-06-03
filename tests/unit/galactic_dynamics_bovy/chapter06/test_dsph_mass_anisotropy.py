"""Unit tests for the dwarf-spheroidal mass-anisotropy degeneracy (Problem 6.6)."""

from pathlib import Path
from unittest.mock import MagicMock, patch

import numpy as np

from galactic_dynamics_bovy.chapter06.dsph_mass_anisotropy import (
    _setup_axes,
    density_normalization,
    enclosed_mass,
    plot_dsph_mass_anisotropy,
    sigma_los_profile,
)

MODULE = "galactic_dynamics_bovy.chapter06.dsph_mass_anisotropy"


class TestEnclosedMass:
    """Tests for the generalized-NFW enclosed mass."""

    def test_increases_with_radius(self) -> None:
        """Enclosed mass grows with radius."""
        m1 = enclosed_mass(0.2, rho_s=6e7, gamma=1.0)
        m2 = enclosed_mass(0.5, rho_s=6e7, gamma=1.0)
        assert 0.0 < m1 < m2

    def test_scales_linearly_with_rho_s(self) -> None:
        """Mass is proportional to rho_s."""
        base = enclosed_mass(0.4, rho_s=1e7, gamma=0.0)
        assert np.isclose(enclosed_mass(0.4, rho_s=3e7, gamma=0.0), 3.0 * base)


class TestDensityNormalization:
    """Tests for matching the enclosed mass across profiles."""

    def test_matches_target_mass(self) -> None:
        """The normalized core encloses the requested mass at r_norm."""
        target = enclosed_mass(0.5, rho_s=6e7, gamma=1.0)
        rho_core = density_normalization(target_mass=target, r_norm=0.5, gamma=0.0)
        assert np.isclose(enclosed_mass(0.5, rho_s=rho_core, gamma=0.0), target, rtol=1e-6)


class TestSigmaLosProfile:
    """Tests for the projected dispersion and the degeneracy."""

    def test_positive_and_finite(self) -> None:
        """The dispersion profile is positive and finite."""
        radii = np.array([0.05, 0.2, 0.5])
        sigma = sigma_los_profile(radii, rho_s=6e7, gamma=1.0, beta=0.0, n_grid=40)
        assert np.all(sigma > 0.0) and np.all(np.isfinite(sigma))

    def test_cored_radial_matches_cusp_better_than_isotropic(self) -> None:
        """A radially anisotropic core fits the cuspy profile far better than an isotropic core."""
        radii = np.array([0.05, 0.1, 0.2, 0.3, 0.45, 0.6])
        target = enclosed_mass(0.5, rho_s=6e7, gamma=1.0)
        rho_core = density_normalization(target_mass=target, r_norm=0.5, gamma=0.0)
        cusp = sigma_los_profile(radii, rho_s=6e7, gamma=1.0, beta=0.0, n_grid=60)
        core_iso = sigma_los_profile(radii, rho_s=rho_core, gamma=0.0, beta=0.0, n_grid=60)
        core_rad = sigma_los_profile(radii, rho_s=rho_core, gamma=0.0, beta=0.5, n_grid=60)
        rms_iso = np.sqrt(np.mean((core_iso - cusp) ** 2))
        rms_rad = np.sqrt(np.mean((core_rad - cusp) ** 2))
        # The radial core is a good fit (well below dSph errors ~2 km/s);
        # the isotropic core is clearly worse.
        assert rms_rad < 0.6
        assert rms_iso > 2.0 * rms_rad

    def test_isotropic_core_rises_outward(self) -> None:
        """An isotropic cored profile predicts a rising sigma_los, unlike the cusp."""
        radii = np.array([0.05, 0.6])
        target = enclosed_mass(0.5, rho_s=6e7, gamma=1.0)
        rho_core = density_normalization(target_mass=target, r_norm=0.5, gamma=0.0)
        sigma = sigma_los_profile(radii, rho_s=rho_core, gamma=0.0, beta=0.0, n_grid=60)
        assert sigma[-1] > sigma[0]


class TestSetupAxes:
    """Tests for the _setup_axes helper."""

    def test_applies_formatting(self) -> None:
        """Should configure ticks."""
        mock_ax = MagicMock()
        _setup_axes(mock_ax)
        assert mock_ax.tick_params.called


class TestPlotDsphMassAnisotropy:
    """Tests for the asset plot function."""

    @patch(f"{MODULE}.sigma_los_profile")
    @patch(f"{MODULE}.plt")
    def test_creates_three_lines(self, mock_plt: MagicMock, mock_profile: MagicMock) -> None:
        """Three sigma_los curves are drawn."""
        mock_ax = MagicMock()
        mock_plt.subplots.return_value = (MagicMock(), mock_ax)
        mock_profile.return_value = np.full(25, 6.5)
        plot_dsph_mass_anisotropy()
        mock_plt.subplots.assert_called_once()
        assert mock_ax.plot.call_count == 3
        assert mock_profile.call_count == 3

    @patch(f"{MODULE}.sigma_los_profile")
    @patch(f"{MODULE}.plt")
    def test_shows_when_no_path(self, mock_plt: MagicMock, mock_profile: MagicMock) -> None:
        """plt.show() is called when no path is given."""
        mock_plt.subplots.return_value = (MagicMock(), MagicMock())
        mock_profile.return_value = np.full(25, 6.5)
        plot_dsph_mass_anisotropy(path=None)
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
        mock_profile.return_value = np.full(25, 6.5)
        plot_dsph_mass_anisotropy(path=tmp_path / "fig.png")
        mock_save.assert_called_once()
        mock_plt.close.assert_called_once_with(mock_fig)
