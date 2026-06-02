"""Unit tests for the NFW Jeans radial dispersion (Problem 5.12 part c)."""

from pathlib import Path
from unittest.mock import MagicMock, patch

import numpy as np

from galactic_dynamics_bovy.chapter05.nfw_jeans_dispersion import (
    _setup_axes,
    dispersion_profile,
    G_KPC_KMS2_MSUN,
    NFWHalo,
    plot_nfw_dispersion,
)

MODULE = "galactic_dynamics_bovy.chapter05.nfw_jeans_dispersion"


class TestNFWHalo:
    """Tests for the NFW halo construction and profiles."""

    def test_from_virial_mass_scales(self) -> None:
        """Milky-Way-like halo has the expected r_vir, r_s and rho_s."""
        halo = NFWHalo.from_virial_mass()
        assert np.isclose(halo.r_vir, 183.2, rtol=1e-2)
        assert np.isclose(halo.r_s, halo.r_vir / 11.5)
        assert halo.rho_s > 0.0

    def test_enclosed_mass_recovers_virial_mass(self) -> None:
        """M(<r_vir) equals the virial mass."""
        halo = NFWHalo.from_virial_mass()
        assert np.isclose(halo.enclosed_mass(halo.r_vir), halo.m_vir, rtol=1e-6)

    def test_density_falls_outward(self) -> None:
        """Density decreases with radius."""
        halo = NFWHalo.from_virial_mass()
        assert halo.density(1.0) > halo.density(10.0) > halo.density(100.0)

    def test_dphi_dr_matches_enclosed_mass(self) -> None:
        """dPhi/dr equals G M(<r) / r^2."""
        halo = NFWHalo.from_virial_mass()
        expected = G_KPC_KMS2_MSUN * halo.enclosed_mass(20.0) / 20.0**2
        assert np.isclose(halo.dphi_dr(20.0), expected)


class TestRadialDispersion:
    """Tests for the direct Jeans integral."""

    def test_isotropic_value(self) -> None:
        """Isotropic dispersion near the scale radius is ~105 km/s."""
        halo = NFWHalo.from_virial_mass()
        sigma = halo.radial_dispersion(halo.r_s, beta=0.0)
        assert np.isclose(sigma, 104.9, rtol=2e-2)

    def test_radial_anisotropy_raises_central_dispersion(self) -> None:
        """At small radii beta=0.5 gives a larger sigma_r than beta=0."""
        halo = NFWHalo.from_virial_mass()
        assert halo.radial_dispersion(1.0, beta=0.5) > halo.radial_dispersion(1.0, beta=0.0)

    def test_dispersion_is_positive(self) -> None:
        """Dispersion is a positive, finite number."""
        halo = NFWHalo.from_virial_mass()
        sigma = halo.radial_dispersion(50.0, beta=0.0)
        assert sigma > 0.0 and np.isfinite(sigma)


class TestDispersionProfile:
    """Tests for the vectorized profile helper."""

    def test_returns_array(self) -> None:
        """Returns one dispersion per input radius."""
        halo = NFWHalo.from_virial_mass()
        radii = np.array([5.0, 20.0, 100.0])
        sigma = dispersion_profile(halo, radii, beta=0.0)
        assert sigma.shape == (3,)
        assert np.all(sigma > 0.0)


class TestSetupAxes:
    """Tests for the _setup_axes helper."""

    def test_applies_formatting(self) -> None:
        """Should configure ticks."""
        mock_ax = MagicMock()
        _setup_axes(mock_ax)
        assert mock_ax.tick_params.called


class TestPlotNFWDispersion:
    """Tests for the asset plot function."""

    @patch(f"{MODULE}.dispersion_profile")
    @patch(f"{MODULE}.plt")
    def test_creates_figure(self, mock_plt: MagicMock, mock_profile: MagicMock) -> None:
        """A figure with both beta curves is created."""
        mock_ax = MagicMock()
        mock_ax.get_ylim.return_value = (0.0, 160.0)
        mock_plt.subplots.return_value = (MagicMock(), mock_ax)
        mock_profile.return_value = np.full(200, 100.0)
        plot_nfw_dispersion()
        mock_plt.subplots.assert_called_once()
        assert mock_ax.semilogx.call_count == 2
        assert mock_profile.call_count == 2

    @patch(f"{MODULE}.dispersion_profile")
    @patch(f"{MODULE}.plt")
    def test_shows_when_no_path(self, mock_plt: MagicMock, mock_profile: MagicMock) -> None:
        """plt.show() is called when no path is given."""
        mock_ax = MagicMock()
        mock_ax.get_ylim.return_value = (0.0, 160.0)
        mock_plt.subplots.return_value = (MagicMock(), mock_ax)
        mock_profile.return_value = np.full(200, 100.0)
        plot_nfw_dispersion(path=None)
        mock_plt.show.assert_called_once()

    @patch(f"{MODULE}.save_figure_if_changed")
    @patch(f"{MODULE}.dispersion_profile")
    @patch(f"{MODULE}.plt")
    def test_saves_when_path_provided(
        self, mock_plt: MagicMock, mock_profile: MagicMock, mock_save: MagicMock, tmp_path: Path
    ) -> None:
        """The figure is saved and closed when a path is provided."""
        mock_fig = MagicMock()
        mock_ax = MagicMock()
        mock_ax.get_ylim.return_value = (0.0, 160.0)
        mock_plt.subplots.return_value = (mock_fig, mock_ax)
        mock_profile.return_value = np.full(200, 100.0)
        plot_nfw_dispersion(path=tmp_path / "fig.png")
        mock_save.assert_called_once()
        mock_plt.close.assert_called_once_with(mock_fig)
