"""Unit tests for the Local Group timing argument (Problem 6.5)."""

from pathlib import Path
from unittest.mock import MagicMock, patch

import numpy as np

from galactic_dynamics_bovy.chapter06.local_group_timing import (
    _setup_axes,
    G_KPC3_MSUN_GYR2,
    KMS_TO_KPC_PER_GYR,
    orbit_track,
    plot_local_group_timing,
    solve_timing_argument,
    TimingArgumentSolution,
)

MODULE = "galactic_dynamics_bovy.chapter06.local_group_timing"


class TestSolveTimingArgument:
    """Tests for the timing-argument solver."""

    def test_local_group_mass(self) -> None:
        """Fiducial inputs give M_LG ~ 4.6e12 Msun."""
        sol = solve_timing_argument()
        assert np.isclose(sol.mass, 4.56e12, rtol=2e-2)

    def test_apocenter(self) -> None:
        """Maximum past separation is ~1050 kpc."""
        sol = solve_timing_argument()
        assert np.isclose(sol.apocenter, 1049.0, rtol=2e-2)
        assert np.isclose(sol.apocenter, 2.0 * sol.semi_major_axis)

    def test_infalling_branch(self) -> None:
        """The eccentric anomaly lies on the infalling branch (pi, 2pi)."""
        sol = solve_timing_argument()
        assert np.pi < sol.eta < 2.0 * np.pi

    def test_reconstructs_velocity(self) -> None:
        """The solution reproduces the input present-day velocity."""
        sol = solve_timing_argument()
        v_kpc_gyr = (
            np.sqrt(G_KPC3_MSUN_GYR2 * sol.mass / sol.semi_major_axis) * np.sin(sol.eta) / (1.0 - np.cos(sol.eta))
        )
        assert np.isclose(v_kpc_gyr / KMS_TO_KPC_PER_GYR, sol.velocity, rtol=1e-6)

    def test_reconstructs_separation(self) -> None:
        """r = a(1 - cos eta) recovers the present separation."""
        sol = solve_timing_argument()
        assert np.isclose(sol.semi_major_axis * (1.0 - np.cos(sol.eta)), sol.separation)

    def test_more_mass_for_faster_approach(self) -> None:
        """A faster approach speed implies a larger inferred mass."""
        slow = solve_timing_argument(velocity=-100.0)
        fast = solve_timing_argument(velocity=-150.0)
        assert fast.mass > slow.mass


class TestOrbitTrack:
    """Tests for the separation history."""

    def test_endpoints(self) -> None:
        """The track runs from r=0 at t=0 to the present separation at the age."""
        sol = solve_timing_argument()
        time, separation = orbit_track(sol, n_points=200)
        assert np.isclose(time[0], 0.0)
        assert np.isclose(separation[0], 0.0)
        assert np.isclose(time[-1], sol.age, rtol=1e-6)
        assert np.isclose(separation[-1], sol.separation, rtol=1e-6)

    def test_apocenter_is_peak(self) -> None:
        """The maximum separation along the track equals the apocenter."""
        sol = solve_timing_argument()
        _, separation = orbit_track(sol, n_points=400)
        assert np.isclose(separation.max(), sol.apocenter, rtol=1e-3)


class TestSetupAxes:
    """Tests for the _setup_axes helper."""

    def test_applies_formatting(self) -> None:
        """Should configure ticks."""
        mock_ax = MagicMock()
        _setup_axes(mock_ax)
        assert mock_ax.tick_params.called


def _mock_solution() -> TimingArgumentSolution:
    """A representative solved timing argument for plot tests."""
    return TimingArgumentSolution(
        eta=4.29,
        semi_major_axis=524.6,
        mass=4.56e12,
        apocenter=1049.0,
        separation=740.0,
        velocity=-125.0,
        age=13.8,
    )


class TestPlotLocalGroupTiming:
    """Tests for the asset plot function."""

    @patch(f"{MODULE}.solve_timing_argument")
    @patch(f"{MODULE}.plt")
    def test_creates_figure(self, mock_plt: MagicMock, mock_solve: MagicMock) -> None:
        """A figure is created from the solved orbit."""
        mock_plt.subplots.return_value = (MagicMock(), MagicMock())
        mock_solve.return_value = _mock_solution()
        plot_local_group_timing()
        mock_plt.subplots.assert_called_once()

    @patch(f"{MODULE}.solve_timing_argument")
    @patch(f"{MODULE}.plt")
    def test_shows_when_no_path(self, mock_plt: MagicMock, mock_solve: MagicMock) -> None:
        """plt.show() is called when no path is given."""
        mock_plt.subplots.return_value = (MagicMock(), MagicMock())
        mock_solve.return_value = _mock_solution()
        plot_local_group_timing(path=None)
        mock_plt.show.assert_called_once()

    @patch(f"{MODULE}.save_figure_if_changed")
    @patch(f"{MODULE}.solve_timing_argument")
    @patch(f"{MODULE}.plt")
    def test_saves_when_path_provided(
        self, mock_plt: MagicMock, mock_solve: MagicMock, mock_save: MagicMock, tmp_path: Path
    ) -> None:
        """The figure is saved and closed when a path is provided."""
        mock_fig = MagicMock()
        mock_plt.subplots.return_value = (mock_fig, MagicMock())
        mock_solve.return_value = _mock_solution()
        plot_local_group_timing(path=tmp_path / "fig.png")
        mock_save.assert_called_once()
        mock_plt.close.assert_called_once_with(mock_fig)
