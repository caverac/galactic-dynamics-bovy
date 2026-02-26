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

| Method | $\Delta\varpi$ (arcmin) |
|---|---|
| 1PN integration (Eqn. 4.11.1) | 12.105 |
| Effective potential quadrature (Problem 4.9) | 12.135 |
| Analytical $6\pi GM_\bullet / [c^2 a(1-e^2)]$ | 12.125 |

All three methods agree to better than 0.3%. The small residual between 1PN integration and the effective-potential result reflects 2PN and higher-order terms: the effective potential derives from the exact Schwarzschild geodesic equation, while the 1PN EOM truncates the post-Newtonian expansion at first order in $GM/(c^2 r)$.

### Orbit comparison

![S2 orbit comparison](assets/generated/p04_11_gr_precession_s2.png)

*Figure 4.11: Radial distance $r(t)$ for S2 over five orbital periods ($T \approx 16.2$ yr). The 1PN orbit (black) progressively leads the Newtonian orbit (gray) due to a shorter radial period. The 1PN correction effectively deepens the potential near pericenter, causing the star to fall in and bounce back faster. Black dots mark 1PN pericenter passages.*

```python
>>> from galactic_dynamics_bovy.chapter04.gr_precession_s2 import plot_gr_precession_s2
>>> plot_gr_precession_s2()
```
