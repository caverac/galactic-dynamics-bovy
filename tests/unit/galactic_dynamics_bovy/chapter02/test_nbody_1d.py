"""Unit tests for 1D N-body simulation module."""

from pathlib import Path
from unittest.mock import MagicMock, patch

import numpy as np

from galactic_dynamics_bovy.chapter02.nbody_1d import (
    compute_forces_1d,
    force_theory_uniform,
    plot_force_comparison,
)


class TestComputeForces1D:
    """Tests for compute_forces_1d function."""

    def test_single_particle(self) -> None:
        """Single particle should feel zero force (no self-interaction)."""
        x = np.array([0.0])
        forces = compute_forces_1d(x, m=1.0, G=1.0)
        # N_left=0, N_right=0, F = 2*pi*(0-0) = 0
        assert np.isclose(forces[0], 0.0)

    def test_two_particles(self) -> None:
        """Two particles should feel opposite forces."""
        x = np.array([-1.0, 1.0])
        forces = compute_forces_1d(x, m=1.0, G=1.0)
        # Particle at -1: N_right=1, N_left=0, F = 2*pi*1*(1-0) = 2*pi
        # Particle at +1: N_right=0, N_left=1, F = 2*pi*1*(0-1) = -2*pi
        assert np.isclose(forces[0], 2 * np.pi)
        assert np.isclose(forces[1], -2 * np.pi)

    def test_three_particles(self) -> None:
        """Middle particle should feel zero net force."""
        x = np.array([-1.0, 0.0, 1.0])
        forces = compute_forces_1d(x, m=1.0, G=1.0)
        # Middle particle: N_right=1, N_left=1, F = 2*pi*(1-1) = 0
        assert np.isclose(forces[1], 0.0)

    def test_preserves_order(self) -> None:
        """Forces should be returned in same order as input positions."""
        x = np.array([1.0, -1.0, 0.0])  # Unsorted
        forces = compute_forces_1d(x, m=1.0, G=1.0)
        # x[0]=1.0 is rightmost, should have negative force
        # x[1]=-1.0 is leftmost, should have positive force
        # x[2]=0.0 is middle, should have zero force
        assert forces[0] < 0
        assert forces[1] > 0
        assert np.isclose(forces[2], 0.0)

    def test_scales_with_mass(self) -> None:
        """Force should scale linearly with mass."""
        x = np.array([-1.0, 1.0])
        f1 = compute_forces_1d(x, m=1.0, G=1.0)
        f2 = compute_forces_1d(x, m=2.0, G=1.0)
        assert np.allclose(f2, 2 * f1)

    def test_scales_with_G(self) -> None:
        """Force should scale linearly with G."""
        x = np.array([-1.0, 1.0])
        f1 = compute_forces_1d(x, m=1.0, G=1.0)
        f2 = compute_forces_1d(x, m=1.0, G=2.0)
        assert np.allclose(f2, 2 * f1)

    def test_sum_of_forces_is_zero(self) -> None:
        """Total force on system should be zero (Newton's third law)."""
        rng = np.random.default_rng(42)
        x = rng.uniform(-1, 1, 100)
        forces = compute_forces_1d(x, m=1.0, G=1.0)
        assert np.isclose(np.sum(forces), 0.0, atol=1e-10)


class TestForceTheoryUniform:
    """Tests for force_theory_uniform function."""

    def test_zero_at_origin(self) -> None:
        """Force should be zero at origin."""
        x = np.array([0.0])
        f = force_theory_uniform(x, total_mass=1.0, G=1.0)
        assert np.isclose(f[0], 0.0)

    def test_antisymmetric(self) -> None:
        """Force should be antisymmetric: F(-x) = -F(x)."""
        x = np.array([0.25])
        f_pos = force_theory_uniform(x, total_mass=1.0, G=1.0)
        f_neg = force_theory_uniform(-x, total_mass=1.0, G=1.0)
        assert np.isclose(f_pos[0], -f_neg[0])

    def test_linear_in_x(self) -> None:
        """Force should be linear in x."""
        x = np.array([0.1, 0.2, 0.3])
        f = force_theory_uniform(x, total_mass=1.0, G=1.0)
        # F = -4*pi*G*A*x, so F/x should be constant
        ratio = f / x
        assert np.allclose(ratio, ratio[0])

    def test_scales_with_mass(self) -> None:
        """Force should scale linearly with total mass."""
        x = np.array([0.25])
        f1 = force_theory_uniform(x, total_mass=1.0, G=1.0)
        f2 = force_theory_uniform(x, total_mass=2.0, G=1.0)
        assert np.isclose(f2[0], 2 * f1[0])

    def test_correct_formula(self) -> None:
        """Check against explicit formula F = -4*pi*G*A*x."""
        x = np.array([0.3])
        A = 5.0
        G = 2.0
        f = force_theory_uniform(x, total_mass=A, G=G)
        expected = -4 * np.pi * G * A * x[0]
        assert np.isclose(f[0], expected)


class TestPlotForceComparison:
    """Tests for plot_force_comparison function."""

    @patch("galactic_dynamics_bovy.chapter02.nbody_1d.plt")
    def test_creates_figure(self, mock_plt: MagicMock) -> None:
        """Verify a figure is created."""
        mock_fig = MagicMock()
        mock_ax = MagicMock()
        mock_plt.subplots.return_value = (mock_fig, mock_ax)

        plot_force_comparison(n_particles=100)

        mock_plt.subplots.assert_called_once()

    @patch("galactic_dynamics_bovy.chapter02.nbody_1d.plt")
    def test_shows_plot_when_no_path(self, mock_plt: MagicMock) -> None:
        """Verify plt.show() is called when no path provided."""
        mock_fig = MagicMock()
        mock_ax = MagicMock()
        mock_plt.subplots.return_value = (mock_fig, mock_ax)

        plot_force_comparison(n_particles=100, path=None)

        mock_plt.show.assert_called_once()

    @patch("galactic_dynamics_bovy.chapter02.nbody_1d.plt")
    def test_saves_figure_when_path_provided(self, mock_plt: MagicMock) -> None:
        """Verify figure is saved when path provided."""
        mock_fig = MagicMock()
        mock_ax = MagicMock()
        mock_plt.subplots.return_value = (mock_fig, mock_ax)
        test_path = Path("/tmp/test_plot.png")

        plot_force_comparison(n_particles=100, path=test_path)

        mock_fig.savefig.assert_called_once()
        mock_plt.close.assert_called_once_with(mock_fig)

    @patch("galactic_dynamics_bovy.chapter02.nbody_1d.plt")
    def test_sets_labels(self, mock_plt: MagicMock) -> None:
        """Verify axis labels are set."""
        mock_fig = MagicMock()
        mock_ax = MagicMock()
        mock_plt.subplots.return_value = (mock_fig, mock_ax)

        plot_force_comparison(n_particles=100)

        mock_ax.set_xlabel.assert_called_once()
        mock_ax.set_ylabel.assert_called_once()


class TestNumericalVsTheory:
    """Integration tests comparing numerical and theoretical forces."""

    def test_ratio_close_to_one(self) -> None:
        """Numerical/theory ratio should be close to 1 for large N."""
        rng = np.random.default_rng(42)
        n = 10000
        x = rng.uniform(-0.5, 0.5, n)

        f_numerical = compute_forces_1d(x, m=1.0, G=1.0)
        f_theory = force_theory_uniform(x, total_mass=n, G=1.0)

        # Exclude points near origin where theory gives small values
        mask = np.abs(x) > 0.1
        ratio = f_numerical[mask] / f_theory[mask]

        # Mean ratio should be close to 1
        assert np.abs(np.mean(ratio) - 1.0) < 0.05
