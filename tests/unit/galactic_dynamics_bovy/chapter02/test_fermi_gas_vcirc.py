"""Unit tests for Fermi gas circular velocity module."""

from pathlib import Path
from unittest.mock import MagicMock, patch

import numpy as np

from galactic_dynamics_bovy.chapter02.fermi_gas_vcirc import (
    get_fermi_gas_params,
    plot_fermi_gas_vcirc,
    vcirc_fermi_gas,
)
from galactic_dynamics_bovy.utils.units import GalacticUnits


class TestGetFermiGasParams:
    """Tests for get_fermi_gas_params function."""

    def test_returns_positive_values(self) -> None:
        """Parameters should be positive."""
        rho_c, k = get_fermi_gas_params(total_mass=10.0, boundary_radius=15.0)
        assert rho_c > 0
        assert k > 0

    def test_k_gives_correct_boundary(self) -> None:
        """k = π/R should give correct boundary radius."""
        boundary_radius = 15.0
        _, k = get_fermi_gas_params(total_mass=10.0, boundary_radius=boundary_radius)
        assert np.isclose(np.pi / k, boundary_radius)

    def test_rho_c_scales_with_mass(self) -> None:
        """Central density should scale linearly with mass."""
        rho_c1, _ = get_fermi_gas_params(total_mass=10.0, boundary_radius=15.0)
        rho_c2, _ = get_fermi_gas_params(total_mass=20.0, boundary_radius=15.0)
        assert np.isclose(rho_c2, 2 * rho_c1)

    def test_rho_c_scales_with_radius_cubed(self) -> None:
        """Central density should scale as R^(-3)."""
        rho_c1, _ = get_fermi_gas_params(total_mass=10.0, boundary_radius=10.0)
        rho_c2, _ = get_fermi_gas_params(total_mass=10.0, boundary_radius=20.0)
        assert np.isclose(rho_c2, rho_c1 / 8)


class TestVcircFermiGas:
    """Tests for vcirc_fermi_gas function."""

    def test_zero_at_origin(self) -> None:
        """Circular velocity should be zero at r=0."""
        units = GalacticUnits()
        rho_c, k = get_fermi_gas_params(10.0, 15.0)
        r = np.array([0.0])
        vc = vcirc_fermi_gas(r, rho_c, k, units.G_kms)
        assert np.isclose(vc[0], 0.0)

    def test_positive_inside_boundary(self) -> None:
        """Circular velocity should be positive inside boundary."""
        units = GalacticUnits()
        rho_c, k = get_fermi_gas_params(10.0, 15.0)
        r = np.array([5.0, 10.0, 14.0])
        vc = vcirc_fermi_gas(r, rho_c, k, units.G_kms)
        assert all(vc > 0)

    def test_finite_at_boundary(self) -> None:
        """Circular velocity should be finite at boundary."""
        units = GalacticUnits()
        boundary_radius = 15.0
        rho_c, k = get_fermi_gas_params(10.0, boundary_radius)
        r = np.array([boundary_radius])
        vc = vcirc_fermi_gas(r, rho_c, k, units.G_kms)
        assert np.isfinite(vc[0])
        assert vc[0] > 0

    def test_scales_with_mass(self) -> None:
        """Velocity should scale as sqrt(M)."""
        units = GalacticUnits()
        r = np.array([10.0])

        rho_c1, k1 = get_fermi_gas_params(10.0, 15.0)
        vc1 = vcirc_fermi_gas(r, rho_c1, k1, units.G_kms)

        rho_c2, k2 = get_fermi_gas_params(40.0, 15.0)
        vc2 = vcirc_fermi_gas(r, rho_c2, k2, units.G_kms)

        assert np.isclose(vc2[0], 2 * vc1[0])

    def test_vectorized(self) -> None:
        """Should work with array input."""
        units = GalacticUnits()
        rho_c, k = get_fermi_gas_params(10.0, 15.0)
        r = np.linspace(0, 15, 100)
        vc = vcirc_fermi_gas(r, rho_c, k, units.G_kms)
        assert len(vc) == 100
        assert all(np.isfinite(vc))


class TestPlotFermiGasVcirc:
    """Tests for plot_fermi_gas_vcirc function."""

    @patch("galactic_dynamics_bovy.chapter02.fermi_gas_vcirc.plt")
    def test_creates_figure(self, mock_plt: MagicMock) -> None:
        """Verify a figure is created."""
        mock_fig = MagicMock()
        mock_ax = MagicMock()
        mock_plt.subplots.return_value = (mock_fig, mock_ax)

        plot_fermi_gas_vcirc()

        mock_plt.subplots.assert_called_once()

    @patch("galactic_dynamics_bovy.chapter02.fermi_gas_vcirc.plt")
    def test_shows_plot_when_no_path(self, mock_plt: MagicMock) -> None:
        """Verify plt.show() is called when no path provided."""
        mock_fig = MagicMock()
        mock_ax = MagicMock()
        mock_plt.subplots.return_value = (mock_fig, mock_ax)

        plot_fermi_gas_vcirc(path=None)

        mock_plt.show.assert_called_once()

    @patch("galactic_dynamics_bovy.chapter02.fermi_gas_vcirc.plt")
    @patch("galactic_dynamics_bovy.chapter02.fermi_gas_vcirc.save_figure_if_changed")
    def test_saves_figure_when_path_provided(
        self,
        mock_save: MagicMock,
        mock_plt: MagicMock,
    ) -> None:
        """Verify figure is saved when path provided."""
        mock_fig = MagicMock()
        mock_ax = MagicMock()
        mock_plt.subplots.return_value = (mock_fig, mock_ax)
        test_path = Path("/tmp/test_plot.png")

        plot_fermi_gas_vcirc(path=test_path)

        mock_save.assert_called_once()
        mock_plt.close.assert_called_once_with(mock_fig)

    @patch("galactic_dynamics_bovy.chapter02.fermi_gas_vcirc.plt")
    def test_sets_labels(self, mock_plt: MagicMock) -> None:
        """Verify axis labels are set."""
        mock_fig = MagicMock()
        mock_ax = MagicMock()
        mock_plt.subplots.return_value = (mock_fig, mock_ax)

        plot_fermi_gas_vcirc()

        mock_ax.set_xlabel.assert_called_once()
        mock_ax.set_ylabel.assert_called_once()
