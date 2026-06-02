"""Unit tests for the cored isothermal sphere profile (Problem 5.12)."""

from pathlib import Path
from unittest.mock import MagicMock, patch

import numpy as np

from galactic_dynamics_bovy.chapter05.cored_isothermal_sphere import (
    _isothermal_rhs,
    _setup_axes,
    cored_isothermal_profile,
    IsothermalProfile,
    plot_cored_isothermal_sphere,
    singular_isothermal_profile,
)

MODULE = "galactic_dynamics_bovy.chapter05.cored_isothermal_sphere"


class TestIsothermalRhs:
    """Tests for the ODE right-hand side."""

    def test_returns_derivatives(self) -> None:
        """At y=0, y'=0 the curvature is y'' = -9 (central value)."""
        deriv = _isothermal_rhs(1.0, np.array([0.0, 0.0]))
        assert deriv[0] == 0.0
        assert np.isclose(deriv[1], -9.0)

    def test_includes_friction_term(self) -> None:
        """The 2/x term enters y''."""
        deriv = _isothermal_rhs(2.0, np.array([0.0, -1.0]))
        # y'' = -9 e^0 - 2*(-1)/2 = -9 + 1 = -8
        assert np.isclose(deriv[1], -8.0)


class TestSingularIsothermalProfile:
    """Tests for the singular isothermal sphere profile."""

    def test_value(self) -> None:
        """f_sing(x) = 2/(9 x^2)."""
        x = np.array([1.0, 2.0])
        result = singular_isothermal_profile(x)
        assert np.allclose(result, [2.0 / 9.0, 2.0 / 36.0])


class TestCoredIsothermalProfile:
    """Tests for the integrated cored isothermal profile (real ODE solve)."""

    def test_returns_profile(self) -> None:
        """Should return an IsothermalProfile with matching array shapes."""
        profile = cored_isothermal_profile(n_points=50)
        assert isinstance(profile, IsothermalProfile)
        assert profile.x.shape == (50,)
        assert profile.f.shape == (50,)

    def test_core_is_flat(self) -> None:
        """Deep in the core f -> 1 (density approaches central value)."""
        profile = cored_isothermal_profile(x_min=1e-2, n_points=100)
        assert np.isclose(profile.f[0], 1.0, atol=1e-3)

    def test_monotonically_decreasing(self) -> None:
        """The density profile decreases outward."""
        profile = cored_isothermal_profile(n_points=200)
        assert np.all(np.diff(profile.f) < 0)

    def test_asymptotes_to_singular(self) -> None:
        """At large x, f x^2 oscillates about the singular value 2/9."""
        profile = cored_isothermal_profile(x_max=1e3, n_points=400)
        tail = profile.f[-1] * profile.x[-1] ** 2
        assert np.isclose(tail, 2.0 / 9.0, atol=0.05)

    def test_core_exceeds_singular(self) -> None:
        """The cored profile is finite where the singular one diverges."""
        profile = cored_isothermal_profile(x_min=1e-2, n_points=50)
        assert profile.f[0] < singular_isothermal_profile(profile.x[0])


class TestSetupAxes:
    """Tests for the _setup_axes helper."""

    def test_applies_formatting(self) -> None:
        """Should configure ticks."""
        mock_ax = MagicMock()
        _setup_axes(mock_ax)
        assert mock_ax.tick_params.called
        mock_ax.grid.assert_not_called()


class TestPlotCoredIsothermalSphere:
    """Tests for the asset plot function."""

    @patch(f"{MODULE}.cored_isothermal_profile")
    @patch(f"{MODULE}.plt")
    def test_creates_figure(self, mock_plt: MagicMock, mock_profile: MagicMock) -> None:
        """A figure with both curves is created."""
        mock_ax = MagicMock()
        mock_plt.subplots.return_value = (MagicMock(), mock_ax)
        x = np.logspace(-2, 2, 20)
        mock_profile.return_value = IsothermalProfile(x=x, f=1.0 / (1.0 + x**2))
        plot_cored_isothermal_sphere()
        mock_plt.subplots.assert_called_once()
        assert mock_ax.loglog.call_count == 2

    @patch(f"{MODULE}.cored_isothermal_profile")
    @patch(f"{MODULE}.plt")
    def test_shows_when_no_path(self, mock_plt: MagicMock, mock_profile: MagicMock) -> None:
        """plt.show() is called when no path is given."""
        mock_plt.subplots.return_value = (MagicMock(), MagicMock())
        x = np.logspace(-2, 2, 20)
        mock_profile.return_value = IsothermalProfile(x=x, f=1.0 / (1.0 + x**2))
        plot_cored_isothermal_sphere(path=None)
        mock_plt.show.assert_called_once()

    @patch(f"{MODULE}.save_figure_if_changed")
    @patch(f"{MODULE}.cored_isothermal_profile")
    @patch(f"{MODULE}.plt")
    def test_saves_when_path_provided(
        self, mock_plt: MagicMock, mock_profile: MagicMock, mock_save: MagicMock, tmp_path: Path
    ) -> None:
        """The figure is saved and closed when a path is provided."""
        mock_fig = MagicMock()
        mock_plt.subplots.return_value = (mock_fig, MagicMock())
        x = np.logspace(-2, 2, 20)
        mock_profile.return_value = IsothermalProfile(x=x, f=1.0 / (1.0 + x**2))
        plot_cored_isothermal_sphere(path=tmp_path / "fig.png")
        mock_save.assert_called_once()
        mock_plt.close.assert_called_once_with(mock_fig)
