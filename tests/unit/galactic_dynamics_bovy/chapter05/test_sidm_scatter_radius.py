"""Unit tests for the SIDM once-scattered radius (Problem 5.12 part d)."""

from pathlib import Path
from unittest.mock import MagicMock, patch

import numpy as np

from galactic_dynamics_bovy.chapter05.nfw_jeans_dispersion import NFWHalo
from galactic_dynamics_bovy.chapter05.sidm_scatter_radius import (
    _setup_axes,
    CROSS_SECTIONS,
    expected_scatters,
    plot_scatter_radius_vs_beta,
    scatter_radius,
    scatter_radius_vs_beta,
)

MODULE = "galactic_dynamics_bovy.chapter05.sidm_scatter_radius"


class TestExpectedScatters:
    """Tests for the expected number of scatters."""

    def test_decreases_outward(self) -> None:
        """Scattering is more frequent in the dense inner halo."""
        halo = NFWHalo.from_virial_mass()
        assert expected_scatters(halo, 5.0, 0.0) > expected_scatters(halo, 50.0, 0.0)

    def test_scales_with_cross_section(self) -> None:
        """Expected scatters is linear in the cross section."""
        halo = NFWHalo.from_virial_mass()
        base = expected_scatters(halo, 10.0, 0.0, cross_section=1.0)
        doubled = expected_scatters(halo, 10.0, 0.0, cross_section=2.0)
        assert np.isclose(doubled, 2.0 * base)


class TestScatterRadius:
    """Tests for the once-scattered radius."""

    def test_isotropic_value(self) -> None:
        """r_1 for beta=0, sigma/m=1, age=10 Gyr is ~16.6 kpc."""
        halo = NFWHalo.from_virial_mass()
        r1 = scatter_radius(halo, 0.0)
        assert np.isclose(r1, 16.6, atol=0.5)

    def test_anisotropic_is_larger(self) -> None:
        """Radial anisotropy raises sigma_r, pushing r_1 outward."""
        halo = NFWHalo.from_virial_mass()
        assert scatter_radius(halo, 0.5) > scatter_radius(halo, 0.0)

    def test_root_satisfies_condition(self) -> None:
        """At r_1 the expected number of scatters is exactly one."""
        halo = NFWHalo.from_virial_mass()
        r1 = scatter_radius(halo, 0.5)
        assert np.isclose(expected_scatters(halo, r1, 0.5), 1.0, atol=1e-6)

    def test_larger_cross_section_grows_radius(self) -> None:
        """A larger cross section thermalizes a larger region."""
        halo = NFWHalo.from_virial_mass()
        assert scatter_radius(halo, 0.0, cross_section=5.0) > scatter_radius(halo, 0.0, cross_section=1.0)


class TestScatterRadiusVsBeta:
    """Tests for the gridded r_1(beta, cross_section) helper."""

    def test_shape_and_consistency(self) -> None:
        """Returns one r_1 per (cross_section, beta) and matches the root finder."""
        halo = NFWHalo.from_virial_mass()
        betas = np.array([0.0, 0.5])
        radii = scatter_radius_vs_beta(halo, betas, (1.0,), n_grid=120)
        assert radii.shape == (1, 2)
        # Interpolated values agree with the direct brentq solution.
        assert np.isclose(radii[0, 0], scatter_radius(halo, 0.0), rtol=2e-2)

    def test_monotonic_in_cross_section(self) -> None:
        """r_1 grows with cross section at fixed beta."""
        halo = NFWHalo.from_virial_mass()
        betas = np.array([0.0])
        radii = scatter_radius_vs_beta(halo, betas, (1.0, 10.0))
        assert radii[1, 0] > radii[0, 0]


class TestSetupAxes:
    """Tests for the _setup_axes helper."""

    def test_applies_formatting(self) -> None:
        """Should configure ticks."""
        mock_ax = MagicMock()
        _setup_axes(mock_ax)
        assert mock_ax.tick_params.called


class TestPlotScatterRadiusVsBeta:
    """Tests for the asset plot function."""

    @patch(f"{MODULE}.scatter_radius_vs_beta")
    @patch(f"{MODULE}.plt")
    def test_creates_figure(self, mock_plt: MagicMock, mock_grid: MagicMock) -> None:
        """A line is drawn for each cross section."""
        mock_ax = MagicMock()
        mock_plt.subplots.return_value = (MagicMock(), mock_ax)
        mock_grid.return_value = np.tile(np.linspace(10.0, 50.0, 30), (len(CROSS_SECTIONS), 1))
        plot_scatter_radius_vs_beta()
        mock_plt.subplots.assert_called_once()
        assert mock_ax.plot.call_count == len(CROSS_SECTIONS)

    @patch(f"{MODULE}.scatter_radius_vs_beta")
    @patch(f"{MODULE}.plt")
    def test_shows_when_no_path(self, mock_plt: MagicMock, mock_grid: MagicMock) -> None:
        """plt.show() is called when no path is given."""
        mock_plt.subplots.return_value = (MagicMock(), MagicMock())
        mock_grid.return_value = np.tile(np.linspace(10.0, 50.0, 30), (len(CROSS_SECTIONS), 1))
        plot_scatter_radius_vs_beta(path=None)
        mock_plt.show.assert_called_once()

    @patch(f"{MODULE}.save_figure_if_changed")
    @patch(f"{MODULE}.scatter_radius_vs_beta")
    @patch(f"{MODULE}.plt")
    def test_saves_when_path_provided(
        self, mock_plt: MagicMock, mock_grid: MagicMock, mock_save: MagicMock, tmp_path: Path
    ) -> None:
        """The figure is saved and closed when a path is provided."""
        mock_fig = MagicMock()
        mock_plt.subplots.return_value = (mock_fig, MagicMock())
        mock_grid.return_value = np.tile(np.linspace(10.0, 50.0, 30), (len(CROSS_SECTIONS), 1))
        plot_scatter_radius_vs_beta(path=tmp_path / "fig.png")
        mock_save.assert_called_once()
        mock_plt.close.assert_called_once_with(mock_fig)
