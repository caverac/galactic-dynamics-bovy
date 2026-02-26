<!-- ======================= -->
<!-- PROBLEM 4.11            -->
<!-- ======================= -->
## Problem 4.11

By directly integrating the first post-Newtonian (1PN) equation of motion

$$
\mathbf{a} = -\frac{GM_\bullet}{r^3}\mathbf{r} + \frac{GM_\bullet}{c^2 r^2}\left(\left[ 4\frac{GM_\bullet}{r^2} - \frac{v^2}{r} \right] \mathbf{r} + 4 v_r \mathbf{v} \right), \tag{4.11.1}
$$

we can compare with the effective-potential quadrature from [Problem 4.9](chapter04.1.md#problem-49) and the analytical approximation. The precession is measured by detecting the next pericenter passage via a zero-crossing event $v_r = \mathbf{r}\cdot\mathbf{v} = 0$ with direction from negative to positive. The precession per orbit is then

$$
\Delta\varpi = \arctan\left(\frac{y_p}{x_p}\right),
$$

where $(x_p, y_p)$ is the position at the next pericenter.

### Three-way comparison

```python
>>> from galactic_dynamics_bovy.chapter04.gr_precession_s2 import compute_precession_1pn, GM_BH, C, A_S2, E_S2
>>> from galactic_dynamics_bovy.chapter04.gr_precession import compute_delta_psi
>>> import numpy as np
>>> prec_1pn = compute_precession_1pn()
>>> r_p = A_S2 * (1.0 - E_S2)
>>> prec_quad = compute_delta_psi(r_p, E_S2, GM=GM_BH, c=C) - 2.0 * np.pi
>>> prec_analytical = 6.0 * np.pi * GM_BH / (C**2 * A_S2 * (1.0 - E_S2**2))
>>> np.degrees(prec_1pn) * 60
12.104969...
>>> np.degrees(prec_quad) * 60
12.135288...
>>> np.degrees(prec_analytical) * 60
12.124627...
```

| Method | $\Delta\varpi$ (arcmin) | vs eff. potential |
|---|---|---|
| Effective potential quadrature ([Problem 4.9](chapter04.1.md#problem-49)) | 12.135 | — |
| 1PN integration (Eqn. 4.11.1) | 12.105 | $-0.25\%$ |
| Analytical $6\pi GM_\bullet / [c^2 a(1-e^2)]$ | 12.125 | $-0.09\%$ |

All three values are consistent with the $\sim 12'$ per orbit measured by the GRAVITY interferometer.

### Why the two GR approaches differ

The effective-potential quadrature from [Problem 4.9](chapter04.1.md#problem-49) adds a correction $-GM_\bullet L^2/(c^2 r^3)$ to the Newtonian effective potential. This term comes from the exact Schwarzschild geodesic equation** — the full metric solution for a test particle around a non-rotating point mass. The corresponding radial force

$$
F(r) = -\frac{GM_\bullet}{r^2} - \frac{3 G M_\bullet L^2}{c^2 r^4}
$$

is conservative and velocity-independent: it can be derived from a scalar potential.

The 1PN equation of motion (Eqn. 4.11.1) comes from a different starting point: it is the test-particle limit of the Einstein-Infeld-Hoffmann (EIH) $N$-body equations [@einstein1938], a weak-field slow-motion expansion ($GM/(c^2 r) \ll 1$, $v/c \ll 1$) of the general-relativistic equations of motion. [@will2008] applies this framework to stellar orbits around Sgr A\*. The EIH formalism is designed for the general $N$-body problem where no exact metric solution is available; Eqn. (4.11.1) is the one-body reduction. The acceleration depends on both position and velocity; in particular the $v_r \mathbf{v}$ term cannot be derived from a scalar potential, reflecting the fact that in GR the gravitational field couples to kinetic energy (mass-energy equivalence).

Both approaches agree at leading order, but because the effective potential captures the exact one-body Schwarzschild result while the 1PN EOM truncates at first order in $GM/(c^2 r)$, the effective-potential quadrature implicitly includes the 2PN and all higher-order corrections that the 1PN integration drops. This is why it gives a slightly larger precession (12.135 vs 12.105 arcmin, a ~0.25% difference). For orbits closer to the black hole the gap would grow.

### Orbit comparison

![S2 orbit comparison](assets/generated/p04_11_gr_precession_s2.png)

*Figure 4.11: Radial distance $r(t)$ for S2 over five orbital periods ($T \approx 16.2$ yr). The 1PN orbit (black) progressively leads the Newtonian orbit (gray) due to a shorter radial period. The 1PN correction effectively deepens the potential near pericenter, causing the star to fall in and bounce back faster. Black dots mark 1PN pericenter passages.*

```python
>>> from galactic_dynamics_bovy.chapter04.gr_precession_s2 import plot_gr_precession_s2
>>> plot_gr_precession_s2()
```
