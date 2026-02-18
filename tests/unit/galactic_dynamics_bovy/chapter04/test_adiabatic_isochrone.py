"""Unit tests for adiabatic isochrone module."""

from pathlib import Path
from unittest.mock import MagicMock, patch

import numpy as np
import numpy.typing as npt

from galactic_dynamics_bovy.chapter04.adiabatic_isochrone import (
    _setup_axes,
    AdiabaticResult,
    estimate_period,
    make_amplitude_func,
    plot_adiabatic_comparison,
    run_adiabatic_experiment,
)


class TestMakeAmplitudeFunc:
    """Tests for make_amplitude_func."""

    def test_returns_callable(self) -> None:
        """Should return a callable."""
        func = make_amplitude_func(tau=1.0, t_mid=5.0)
        assert callable(func)

    def test_asymptotic_low(self) -> None:
        """For t << t_mid, amplitude should be close to 1."""
        func = make_amplitude_func(tau=1.0, t_mid=50.0)
        assert np.isclose(func(0.0), 1.0, atol=0.01)

    def test_asymptotic_high(self) -> None:
        """For t >> t_mid, amplitude should be close to 2."""
        func = make_amplitude_func(tau=1.0, t_mid=50.0)
        assert np.isclose(func(100.0), 2.0, atol=0.01)

    def test_midpoint_value(self) -> None:
        """At t = t_mid, amplitude should be 1.5."""
        func = make_amplitude_func(tau=1.0, t_mid=50.0)
        assert np.isclose(func(50.0), 1.5, atol=1e-10)

    def test_monotonically_increasing(self) -> None:
        """Amplitude should be monotonically increasing."""
        func = make_amplitude_func(tau=1.0, t_mid=50.0)
        t_vals = np.linspace(0.0, 100.0, 100)
        a_vals = [func(t) for t in t_vals]
        diffs = np.diff(a_vals)
        assert np.all(diffs >= 0)

    def test_returns_float(self) -> None:
        """Should return a float."""
        func = make_amplitude_func(tau=1.0, t_mid=5.0)
        result = func(3.0)
        assert isinstance(result, float)


class TestEstimatePeriod:
    """Tests for estimate_period."""

    def test_returns_positive(self) -> None:
        """Period should be positive."""
        T = estimate_period(b=1.0)
        assert T > 0

    def test_reasonable_magnitude(self) -> None:
        """Period should be order-of-magnitude reasonable for natural units."""
        T = estimate_period(b=1.0)
        assert 1.0 < T < 100.0


class TestRunAdiabaticExperiment:
    """Tests for run_adiabatic_experiment."""

    def test_returns_adiabatic_result(self) -> None:
        """Should return an AdiabaticResult dataclass."""
        result = run_adiabatic_experiment(tau=1.0, n_orbits=5, n_steps=1001, n_snapshots=20)
        assert isinstance(result, AdiabaticResult)

    def test_output_shapes(self) -> None:
        """Output arrays should have consistent shapes."""
        n_steps = 1001
        n_snapshots = 20
        result = run_adiabatic_experiment(tau=1.0, n_orbits=5, n_steps=n_steps, n_snapshots=n_snapshots)
        assert len(result.t) == n_steps
        assert len(result.R) == n_steps
        assert len(result.Jr) == n_snapshots
        assert len(result.Lz) == n_snapshots
        assert len(result.amplitude) == n_snapshots

    def test_amplitude_ramps(self) -> None:
        """Amplitude should start near 1 and end near 2."""
        result = run_adiabatic_experiment(tau=0.5, n_orbits=10, n_steps=1001, n_snapshots=50)
        assert np.isclose(result.amplitude[0], 1.0, atol=0.05)
        assert np.isclose(result.amplitude[-1], 2.0, atol=0.05)

    def test_positive_Jr(self) -> None:
        """Radial action should be positive (eccentric orbit)."""
        result = run_adiabatic_experiment(tau=1.0, n_orbits=5, n_steps=1001, n_snapshots=20)
        assert np.all(result.Jr > 0)

    def test_tau_stored(self) -> None:
        """Tau should be stored in the result."""
        result = run_adiabatic_experiment(tau=3.14, n_orbits=5, n_steps=1001, n_snapshots=20)
        assert result.tau == 3.14

    def test_T_orbit_stored(self) -> None:
        """T_orbit should be stored and positive."""
        result = run_adiabatic_experiment(tau=1.0, n_orbits=5, n_steps=1001, n_snapshots=20)
        assert result.T_orbit > 0

    def test_slow_ramp_preserves_Jr(self) -> None:
        """Slow ramp should approximately preserve Jr (adiabatic invariant)."""
        T_orbit = estimate_period()
        tau_slow = 10.0 * T_orbit
        result = run_adiabatic_experiment(tau_slow, n_orbits=60, n_steps=20001, n_snapshots=200)
        # Compare Jr before and after the ramp
        # Use early and late values (avoid edges where ramp is happening)
        jr_early = np.median(result.Jr[:20])
        jr_late = np.median(result.Jr[-20:])
        fractional_change = abs(jr_late - jr_early) / jr_early
        assert fractional_change < 0.15

    def test_fast_ramp_changes_Jr(self) -> None:
        """Fast ramp should change Jr significantly."""
        T_orbit = estimate_period()
        tau_fast = 0.1 * T_orbit
        result = run_adiabatic_experiment(tau_fast, n_orbits=40, n_steps=10001, n_snapshots=200)
        jr_early = np.median(result.Jr[:20])
        jr_late = np.median(result.Jr[-20:])
        fractional_change = abs(jr_late - jr_early) / jr_early
        assert fractional_change > 0.15


class TestSetupAxes:
    """Tests for _setup_axes helper."""

    def test_applies_formatting(self) -> None:
        """Should call tick_params and grid on the axes."""
        mock_ax = MagicMock()
        _setup_axes(mock_ax)
        assert mock_ax.tick_params.called
        mock_ax.grid.assert_called_once()


def _make_mock_result() -> AdiabaticResult:
    """Create a mock AdiabaticResult for plot tests."""
    n_steps = 100
    n_snap = 50
    return AdiabaticResult(
        t=np.linspace(0, 100, n_steps),
        R=np.ones(n_steps),
        Jr=np.full(n_snap, 0.1),
        Lz=np.full(n_snap, 0.5),
        amplitude=np.linspace(1.0, 2.0, n_snap),
        T_orbit=5.0,
        tau=1.0,
    )


def _make_mock_2x2_axes() -> tuple[npt.NDArray[np.object_], list[MagicMock]]:
    """Create a 2x2 numpy object array of MagicMock axes.

    Returns the array and a flat list of the 4 mock axes for assertions.
    """
    axes = np.empty((2, 2), dtype=object)
    flat: list[MagicMock] = []
    for i in range(2):
        for j in range(2):
            m = MagicMock()
            m.get_ylim.return_value = (0.0, 1.0)
            axes[i, j] = m
            flat.append(m)
    return axes, flat


class TestPlotAdiabaticComparison:
    """Tests for plot_adiabatic_comparison function."""

    @patch("galactic_dynamics_bovy.chapter04.adiabatic_isochrone.run_adiabatic_experiment")
    @patch("galactic_dynamics_bovy.chapter04.adiabatic_isochrone.estimate_period")
    @patch("galactic_dynamics_bovy.chapter04.adiabatic_isochrone.plt")
    def test_creates_figure(self, mock_plt: MagicMock, mock_period: MagicMock, mock_run: MagicMock) -> None:
        """Verify a 2x2 figure is created."""
        mock_fig = MagicMock()
        mock_axes, _ = _make_mock_2x2_axes()
        mock_plt.subplots.return_value = (mock_fig, mock_axes)
        mock_period.return_value = 5.0
        mock_run.return_value = _make_mock_result()
        plot_adiabatic_comparison()
        mock_plt.subplots.assert_called_once()

    @patch("galactic_dynamics_bovy.chapter04.adiabatic_isochrone.run_adiabatic_experiment")
    @patch("galactic_dynamics_bovy.chapter04.adiabatic_isochrone.estimate_period")
    @patch("galactic_dynamics_bovy.chapter04.adiabatic_isochrone.plt")
    def test_shows_plot_when_no_path(self, mock_plt: MagicMock, mock_period: MagicMock, mock_run: MagicMock) -> None:
        """Verify plt.show() is called when no path provided."""
        mock_axes, _ = _make_mock_2x2_axes()
        mock_plt.subplots.return_value = (MagicMock(), mock_axes)
        mock_period.return_value = 5.0
        mock_run.return_value = _make_mock_result()
        plot_adiabatic_comparison(path=None)
        mock_plt.show.assert_called_once()

    @patch("galactic_dynamics_bovy.chapter04.adiabatic_isochrone.run_adiabatic_experiment")
    @patch("galactic_dynamics_bovy.chapter04.adiabatic_isochrone.estimate_period")
    @patch("galactic_dynamics_bovy.chapter04.adiabatic_isochrone.plt")
    def test_saves_figure_when_path_provided(
        self, mock_plt: MagicMock, mock_period: MagicMock, mock_run: MagicMock
    ) -> None:
        """Verify figure is saved when path provided."""
        mock_fig = MagicMock()
        mock_axes, _ = _make_mock_2x2_axes()
        mock_plt.subplots.return_value = (mock_fig, mock_axes)
        mock_period.return_value = 5.0
        mock_run.return_value = _make_mock_result()
        test_path = Path("/tmp/test_adiabatic.png")
        plot_adiabatic_comparison(path=test_path)
        mock_fig.savefig.assert_called_once()
        mock_plt.close.assert_called_once_with(mock_fig)

    @patch("galactic_dynamics_bovy.chapter04.adiabatic_isochrone.run_adiabatic_experiment")
    @patch("galactic_dynamics_bovy.chapter04.adiabatic_isochrone.estimate_period")
    @patch("galactic_dynamics_bovy.chapter04.adiabatic_isochrone.plt")
    def test_sets_labels(self, mock_plt: MagicMock, mock_period: MagicMock, mock_run: MagicMock) -> None:
        """Verify xlabel on bottom row, ylabel on left column only."""
        mock_axes, flat = _make_mock_2x2_axes()
        mock_plt.subplots.return_value = (MagicMock(), mock_axes)
        mock_period.return_value = 5.0
        mock_run.return_value = _make_mock_result()
        plot_adiabatic_comparison()
        # Bottom row gets xlabel
        flat[2].set_xlabel.assert_called_once()
        flat[3].set_xlabel.assert_called_once()
        # Top row does not
        flat[0].set_xlabel.assert_not_called()
        flat[1].set_xlabel.assert_not_called()
        # Left column gets ylabel
        flat[0].set_ylabel.assert_called_once()
        flat[2].set_ylabel.assert_called_once()
        # Right column does not
        flat[1].set_ylabel.assert_not_called()
        flat[3].set_ylabel.assert_not_called()

    @patch("galactic_dynamics_bovy.chapter04.adiabatic_isochrone.run_adiabatic_experiment")
    @patch("galactic_dynamics_bovy.chapter04.adiabatic_isochrone.estimate_period")
    @patch("galactic_dynamics_bovy.chapter04.adiabatic_isochrone.plt")
    def test_runs_two_experiments(self, mock_plt: MagicMock, mock_period: MagicMock, mock_run: MagicMock) -> None:
        """Should call run_adiabatic_experiment twice (fast and slow)."""
        mock_axes, _ = _make_mock_2x2_axes()
        mock_plt.subplots.return_value = (MagicMock(), mock_axes)
        mock_period.return_value = 5.0
        mock_run.return_value = _make_mock_result()
        plot_adiabatic_comparison()
        assert mock_run.call_count == 2
