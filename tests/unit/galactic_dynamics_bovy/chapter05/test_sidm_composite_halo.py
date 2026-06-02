"""Unit tests for the composite SIDM halo model (Problem 5.12 part e)."""

from pathlib import Path
from unittest.mock import MagicMock, patch

import numpy as np

from galactic_dynamics_bovy.chapter05.nfw_jeans_dispersion import NFWHalo
from galactic_dynamics_bovy.chapter05.sidm_composite_halo import (
    _setup_axes,
    composite_density,
    CoredMatch,
    match_cored_isothermal_to_nfw,
    plot_sidm_composite_halo,
)

MODULE = "galactic_dynamics_bovy.chapter05.sidm_composite_halo"

# beta = 0 once-scattered radius for the default halo (Problem 5.12d).
R1 = 16.576


class TestMatchCoredIsothermalToNFW:
    """Tests for the density-and-mass matching at r_1."""

    def test_fitted_parameters(self) -> None:
        """The fitted rho_0, r_0 and X_1 reproduce the prototype values."""
        halo = NFWHalo.from_virial_mass()
        match = match_cored_isothermal_to_nfw(halo, R1)
        assert np.isclose(match.x_match, 1.995, rtol=2e-2)
        assert np.isclose(match.r_0, 8.31, rtol=2e-2)
        assert np.isclose(match.rho_0, 2.32e7, rtol=3e-2)

    def test_density_matches_nfw_at_r1(self) -> None:
        """Cored density equals NFW density at the matching radius."""
        halo = NFWHalo.from_virial_mass()
        match = match_cored_isothermal_to_nfw(halo, R1)
        rho_cored = composite_density(halo, match, np.array([R1 * 0.999]))[0]
        assert np.isclose(rho_cored, halo.density(R1), rtol=3e-3)

    def test_core_is_denser_than_center_value(self) -> None:
        """The cored central density exceeds the NFW value at r_1 (flat core)."""
        halo = NFWHalo.from_virial_mass()
        match = match_cored_isothermal_to_nfw(halo, R1)
        assert match.rho_0 > halo.density(R1)


class TestCompositeDensity:
    """Tests for the piecewise composite density."""

    def test_uses_nfw_outside_r1(self) -> None:
        """Outside r_1 the composite equals the NFW density."""
        halo = NFWHalo.from_virial_mass()
        match = match_cored_isothermal_to_nfw(halo, R1)
        radii = np.array([50.0, 100.0])
        composite = composite_density(halo, match, radii)
        nfw = np.array([halo.density(float(r)) for r in radii])
        assert np.allclose(composite, nfw)

    def test_continuous_at_r1(self) -> None:
        """The matched cored density equals the NFW density at r_1."""
        halo = NFWHalo.from_virial_mass()
        match = match_cored_isothermal_to_nfw(halo, R1)
        value = composite_density(halo, match, np.array([R1]))[0]
        assert np.isclose(value, halo.density(R1), rtol=1e-6)

    def test_core_flattens_inward(self) -> None:
        """Inside r_1 the cored density stays below the NFW cusp near the centre."""
        halo = NFWHalo.from_virial_mass()
        match = match_cored_isothermal_to_nfw(halo, R1)
        composite = composite_density(halo, match, np.array([1.0]))[0]
        assert composite < halo.density(1.0)


class TestSetupAxes:
    """Tests for the _setup_axes helper."""

    def test_applies_formatting(self) -> None:
        """Should configure ticks."""
        mock_ax = MagicMock()
        _setup_axes(mock_ax)
        assert mock_ax.tick_params.called


class TestPlotSidmCompositeHalo:
    """Tests for the asset plot function."""

    @patch(f"{MODULE}.composite_density")
    @patch(f"{MODULE}.match_cored_isothermal_to_nfw")
    @patch(f"{MODULE}.scatter_radius")
    @patch(f"{MODULE}.plt")
    def test_creates_figure(
        self, mock_plt: MagicMock, mock_r1: MagicMock, mock_match: MagicMock, mock_density: MagicMock
    ) -> None:
        """Both the NFW and composite curves are drawn."""
        mock_ax = MagicMock()
        mock_ax.get_ylim.return_value = (1e3, 1e8)
        mock_plt.subplots.return_value = (MagicMock(), mock_ax)
        mock_r1.return_value = R1
        mock_match.return_value = CoredMatch(rho_0=2.3e7, r_0=8.3, r_1=R1, x_match=2.0)
        mock_density.return_value = np.full(300, 1e6)
        plot_sidm_composite_halo()
        mock_plt.subplots.assert_called_once()
        assert mock_ax.loglog.call_count == 2

    @patch(f"{MODULE}.composite_density")
    @patch(f"{MODULE}.match_cored_isothermal_to_nfw")
    @patch(f"{MODULE}.scatter_radius")
    @patch(f"{MODULE}.plt")
    def test_shows_when_no_path(
        self, mock_plt: MagicMock, mock_r1: MagicMock, mock_match: MagicMock, mock_density: MagicMock
    ) -> None:
        """plt.show() is called when no path is given."""
        mock_ax = MagicMock()
        mock_ax.get_ylim.return_value = (1e3, 1e8)
        mock_plt.subplots.return_value = (MagicMock(), mock_ax)
        mock_r1.return_value = R1
        mock_match.return_value = CoredMatch(rho_0=2.3e7, r_0=8.3, r_1=R1, x_match=2.0)
        mock_density.return_value = np.full(300, 1e6)
        plot_sidm_composite_halo(path=None)
        mock_plt.show.assert_called_once()

    @patch(f"{MODULE}.save_figure_if_changed")
    @patch(f"{MODULE}.composite_density")
    @patch(f"{MODULE}.match_cored_isothermal_to_nfw")
    @patch(f"{MODULE}.scatter_radius")
    @patch(f"{MODULE}.plt")
    def test_saves_when_path_provided(  # pylint: disable=too-many-positional-arguments
        self,
        mock_plt: MagicMock,
        mock_r1: MagicMock,
        mock_match: MagicMock,
        mock_density: MagicMock,
        mock_save: MagicMock,
        tmp_path: Path,
    ) -> None:
        """The figure is saved and closed when a path is provided."""
        mock_fig = MagicMock()
        mock_ax = MagicMock()
        mock_ax.get_ylim.return_value = (1e3, 1e8)
        mock_plt.subplots.return_value = (mock_fig, mock_ax)
        mock_r1.return_value = R1
        mock_match.return_value = CoredMatch(rho_0=2.3e7, r_0=8.3, r_1=R1, x_match=2.0)
        mock_density.return_value = np.full(300, 1e6)
        plot_sidm_composite_halo(path=tmp_path / "fig.png")
        mock_save.assert_called_once()
        mock_plt.close.assert_called_once_with(mock_fig)
