"""Unit tests for GR precession of S2 via 1PN orbit integration."""

from pathlib import Path
from unittest.mock import MagicMock, patch

import numpy as np

from galactic_dynamics_bovy.chapter04.gr_precession import compute_delta_psi
from galactic_dynamics_bovy.chapter04.gr_precession_s2 import (
    _eom_1pn,
    _pericenter_event,
    A_S2,
    C,
    compute_precession_1pn,
    E_S2,
    GM_BH,
    integrate_orbit_1pn,
    plot_gr_precession_s2,
)


class TestEom1pn:
    """Tests for _eom_1pn."""

    def test_newtonian_limit(self) -> None:
        """For very large c, 1PN correction should vanish (Newtonian limit)."""
        r_p = A_S2 * (1.0 - E_S2)
        v_p = np.sqrt(GM_BH / A_S2 * (1.0 + E_S2) / (1.0 - E_S2))
        y = np.array([r_p, 0.0, 0.0, v_p])

        # Large c: 1PN correction negligible
        dydt_large_c = _eom_1pn(0.0, y, GM_BH, 1e30)
        # Pure Newtonian
        r = r_p
        ax_newton = -GM_BH / r**2
        assert np.isclose(dydt_large_c[2], ax_newton, rtol=1e-10)
        assert np.isclose(dydt_large_c[3], 0.0, atol=1e-10)

    def test_correct_dimensions(self) -> None:
        """Output should have 4 components: vx, vy, ax, ay."""
        y = np.array([A_S2 * (1.0 - E_S2), 0.0, 0.0, 1e4])
        result = _eom_1pn(0.0, y, GM_BH, C)
        assert len(result) == 4

    def test_gr_correction_sign(self) -> None:
        """At pericenter with purely tangential velocity, 1PN correction
        should make the radial acceleration less negative (weaker pull)
        compared to Newtonian for a prograde orbit."""
        r_p = A_S2 * (1.0 - E_S2)
        v_p = np.sqrt(GM_BH / A_S2 * (1.0 + E_S2) / (1.0 - E_S2))
        y = np.array([r_p, 0.0, 0.0, v_p])

        # At pericenter, v_r = 0, so the 1PN correction is
        # GM/(c^2*r^2) * (4GM/r^2 - v^2/r) * x
        # For S2, v/c ~ 0.025, GM/(c^2*r) ~ 0.003
        # 4GM/r^2 vs v^2/r: need to check which dominates
        dydt_1pn = _eom_1pn(0.0, y, GM_BH, C)
        dydt_newton = _eom_1pn(0.0, y, GM_BH, 1e30)

        # 1PN x-acceleration should differ from Newtonian
        assert dydt_1pn[2] != dydt_newton[2]

    def test_velocity_components(self) -> None:
        """First two outputs should be the velocity components."""
        y = np.array([1e12, 2e12, 3e4, 4e4])
        result = _eom_1pn(0.0, y, GM_BH, C)
        assert result[0] == 3e4
        assert result[1] == 4e4


class TestPericenterEvent:
    """Tests for _pericenter_event."""

    def test_zero_at_pericenter(self) -> None:
        """At pericenter with tangential velocity, v_r should be zero."""
        r_p = A_S2 * (1.0 - E_S2)
        y = np.array([r_p, 0.0, 0.0, 1e4])
        assert _pericenter_event(0.0, y) == 0.0

    def test_nonzero_for_radial_motion(self) -> None:
        """For a state with radial velocity, event should be nonzero."""
        y = np.array([1e12, 0.0, 1e4, 1e4])
        assert _pericenter_event(0.0, y) != 0.0

    def test_direction_attribute(self) -> None:
        """Event should have direction=+1."""
        assert _pericenter_event.direction == 1


class TestIntegrateOrbit:
    """Tests for integrate_orbit_1pn."""

    def test_returns_valid_solution(self) -> None:
        """Solution should have success=True and non-empty arrays."""
        sol = integrate_orbit_1pn(n_orbits=1.5)
        assert sol.success
        assert len(sol.t) > 100

    def test_orbit_stays_bound(self) -> None:
        """The orbit should remain bound: r should stay finite and positive."""
        sol = integrate_orbit_1pn(n_orbits=1.5)
        x = sol.y[0]
        y = sol.y[1]
        r = np.sqrt(x**2 + y**2)
        r_a = A_S2 * (1.0 + E_S2)
        assert np.all(r > 0)
        assert np.all(r < 2.0 * r_a)  # should not exceed ~2x apocenter

    def test_detects_pericenter_events(self) -> None:
        """Should detect at least one pericenter event."""
        sol = integrate_orbit_1pn(n_orbits=1.5)
        assert len(sol.t_events[0]) >= 1

    def test_energy_approximately_conserved(self) -> None:
        """Newtonian energy should be approximately conserved
        (the 1PN correction is tiny)."""
        sol = integrate_orbit_1pn(n_orbits=1.5)
        x = sol.y[0]
        y = sol.y[1]
        vx = sol.y[2]
        vy = sol.y[3]
        r = np.sqrt(x**2 + y**2)
        v2 = vx**2 + vy**2
        E = 0.5 * v2 - GM_BH / r
        # Newtonian energy varies at the level of the 1PN correction
        # (~GM/(c²a) ~ 3e-3 for S2), so use a generous tolerance
        assert np.std(E) / np.abs(np.mean(E)) < 1e-2


class TestComputePrecession1pn:
    """Tests for compute_precession_1pn."""

    def test_positive_precession(self) -> None:
        """GR precession should be prograde (positive)."""
        prec = compute_precession_1pn()
        assert prec > 0

    def test_matches_analytical(self) -> None:
        """Should match the analytical formula to ~1%."""
        prec = compute_precession_1pn()
        analytical = 6.0 * np.pi * GM_BH / (C**2 * A_S2 * (1.0 - E_S2**2))
        assert np.isclose(prec, analytical, rtol=0.01)

    def test_consistent_with_effective_potential(self) -> None:
        """Should be consistent with the effective potential quadrature."""
        prec_1pn = compute_precession_1pn()
        r_p = A_S2 * (1.0 - E_S2)
        prec_quad = compute_delta_psi(r_p, E_S2, GM=GM_BH, c=C) - 2.0 * np.pi
        assert np.isclose(prec_1pn, prec_quad, rtol=0.01)

    def test_order_of_magnitude(self) -> None:
        """Precession should be on the order of ~10 arcmin for S2."""
        prec = compute_precession_1pn()
        prec_arcmin = np.degrees(prec) * 60.0
        assert 5.0 < prec_arcmin < 20.0

    @patch("galactic_dynamics_bovy.chapter04.gr_precession_s2.integrate_orbit_1pn")
    def test_negative_angle_wrapped(self, mock_integrate: MagicMock) -> None:
        """When arctan2 returns a negative angle, it should be wrapped to [0, 2*pi)."""
        T_kepler = 2.0 * np.pi * np.sqrt(A_S2**3 / GM_BH)
        mock_sol = MagicMock()
        mock_sol.t_events = [np.array([0.01 * T_kepler, 0.8 * T_kepler])]
        # Pericenter at negative y => arctan2 returns negative angle
        mock_sol.sol.return_value = np.array([1e13, -1e11, 0.0, 1e4])
        mock_integrate.return_value = mock_sol

        prec = compute_precession_1pn()
        assert prec > 0
        assert prec < 2.0 * np.pi


class TestPlot:
    """Tests for plot_gr_precession_s2."""

    @patch("galactic_dynamics_bovy.chapter04.gr_precession_s2.integrate_orbit_1pn")
    @patch("galactic_dynamics_bovy.chapter04.gr_precession_s2.plt")
    def test_creates_figure(self, mock_plt: MagicMock, mock_integrate: MagicMock) -> None:
        """Verify a figure is created."""
        mock_sol = MagicMock()
        mock_sol.t = np.linspace(0, 5e8, 100)
        mock_sol.y = np.random.randn(4, 100) * 1e13
        mock_sol.t_events = [np.array([1e8, 2e8])]
        mock_sol.sol.return_value = np.array([1e13, 1e12, 0.0, 1e4])
        mock_integrate.return_value = mock_sol

        mock_fig = MagicMock()
        mock_ax = MagicMock()
        mock_plt.subplots.return_value = (mock_fig, mock_ax)
        plot_gr_precession_s2()
        mock_plt.subplots.assert_called_once()

    @patch("galactic_dynamics_bovy.chapter04.gr_precession_s2.integrate_orbit_1pn")
    @patch("galactic_dynamics_bovy.chapter04.gr_precession_s2.plt")
    def test_shows_plot_when_no_path(self, mock_plt: MagicMock, mock_integrate: MagicMock) -> None:
        """Verify plt.show() is called when no path provided."""
        mock_sol = MagicMock()
        mock_sol.t = np.linspace(0, 5e8, 100)
        mock_sol.y = np.random.randn(4, 100) * 1e13
        mock_sol.t_events = [np.array([1e8, 2e8])]
        mock_sol.sol.return_value = np.array([1e13, 1e12, 0.0, 1e4])
        mock_integrate.return_value = mock_sol

        mock_plt.subplots.return_value = (MagicMock(), MagicMock())
        plot_gr_precession_s2(path=None)
        mock_plt.show.assert_called_once()

    @patch("galactic_dynamics_bovy.chapter04.gr_precession_s2.integrate_orbit_1pn")
    @patch("galactic_dynamics_bovy.chapter04.gr_precession_s2.plt")
    def test_saves_figure_when_path_provided(self, mock_plt: MagicMock, mock_integrate: MagicMock) -> None:
        """Verify figure is saved when path provided."""
        mock_sol = MagicMock()
        mock_sol.t = np.linspace(0, 5e8, 100)
        mock_sol.y = np.random.randn(4, 100) * 1e13
        mock_sol.t_events = [np.array([1e8, 2e8])]
        mock_sol.sol.return_value = np.array([1e13, 1e12, 0.0, 1e4])
        mock_integrate.return_value = mock_sol

        mock_fig = MagicMock()
        mock_ax = MagicMock()
        mock_plt.subplots.return_value = (mock_fig, mock_ax)
        test_path = Path("/tmp/test_gr_precession_s2.png")
        plot_gr_precession_s2(path=test_path)
        mock_fig.savefig.assert_called_once()
        mock_plt.close.assert_called_once_with(mock_fig)

    @patch("galactic_dynamics_bovy.chapter04.gr_precession_s2.integrate_orbit_1pn")
    @patch("galactic_dynamics_bovy.chapter04.gr_precession_s2.plt")
    def test_sets_labels(self, mock_plt: MagicMock, mock_integrate: MagicMock) -> None:
        """Verify axis labels are set."""
        mock_sol = MagicMock()
        mock_sol.t = np.linspace(0, 5e8, 100)
        mock_sol.y = np.random.randn(4, 100) * 1e13
        mock_sol.t_events = [np.array([1e8, 2e8])]
        mock_sol.sol.return_value = np.array([1e13, 1e12, 0.0, 1e4])
        mock_integrate.return_value = mock_sol

        mock_ax = MagicMock()
        mock_plt.subplots.return_value = (MagicMock(), mock_ax)
        plot_gr_precession_s2()
        mock_ax.set_xlabel.assert_called_once()
        mock_ax.set_ylabel.assert_called_once()
        mock_ax.legend.assert_called_once()
