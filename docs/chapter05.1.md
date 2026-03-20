# Chapter 5: Equilibria of collisionless stellar systems

<!-- ======================= -->
<!-- PROBLEM 5.1             -->
<!-- ======================= -->
## Problem 5.1

Assuming each star is roughly $1\,M_\odot$, Fornax has $N = 10^7$ stars. The relaxation time is approximately

$$
\begin{align}
t_{\mathrm{relax}} &\approx \frac{N}{8\ln N}\,t_{\mathrm{cross}} \\
&= \frac{N}{8\ln N}\,\frac{2\pi R}{\sigma} \\
&\approx \frac{10^7}{8 \times 16.1} \times \frac{2\pi \times 1\,\mathrm{kpc}}{10\,\mathrm{km/s}} \\
&\approx 4.7 \times 10^4\,\mathrm{Gyr}
\end{align}
$$

<!-- ======================= -->
<!-- PROBLEM 5.2             -->
<!-- ======================= -->
## Problem 5.2

From the virial theorem $2K + W = 0$, if the mass increases by a factor of $f$ and all the distances remain the same, then $W$ increases by a factor of $f^2$ (since $W \propto M^2/R$) and thus $K$ must also increase by a factor of $f^2$. Since $K \propto M \sigma^2$, $\sigma$ must increase by a factor of $\sqrt{f}$.

<!-- ======================= -->
<!-- PROBLEM 5.3             -->
<!-- ======================= -->
## Problem 5.3

### Acceleration is dominated by distant bodies

The gravitational force from a single body at distance $r$ is $Gm/r^2$. In a roughly homogeneous system the enclosed mass grows as $M(<r) \propto r^3$, so the cumulative acceleration from all mass inside radius $r$ is

$$
a(<r) \sim \frac{GM(<r)}{r^2} \propto r
$$

which increases with $r$: distant matter dominates. More concretely, the nearest neighbor sits at $d \sim R/N^{1/3}$ and contributes $a_{\mathrm{near}} \sim Gm/d^2$, while the total smooth-field acceleration is $a_{\mathrm{smooth}} \sim GNm/R^2$. Their ratio is

$$
\frac{a_{\mathrm{near}}}{a_{\mathrm{smooth}}} \sim \frac{1}{N}\left(\frac{R}{d}\right)^2 \sim N^{-1/3} \ll 1
$$

so the nearest neighbor's contribution is negligible.

### Jerk is dominated by nearby bodies

The jerk from a body at distance $r$ involves the tidal tensor $\partial^2\Phi/\partial x_i\partial x_j \sim Gm/r^3$, giving a jerk that scales as $Gmv/r^3$. The contribution from particles in a shell of radius $r$ and thickness $dr$ to the total stochastic jerk magnitude is

$$
d\dot{a}(r) \sim n\,r^2 dr \cdot \frac{Gmv}{r^3} = Gnmv \frac{dr}{r}
$$

which implies $\dot{a}$ logarithmically diverges as $r \to 0$: nearby bodies dominate. The nearest neighbor at distance $d \sim R/N^{1/3}$ contributes

$$
\dot{a}_{\mathrm{near}} \sim \frac{Gmv}{d^3} \sim Gnmv \sim G\bar{\rho}\,v
$$

which is of the same order as the jerk from the entire smooth-field tidal tensor ($\dot{a}_{\mathrm{smooth}} \sim G\bar{\rho}\,v$). Unlike the acceleration, where the nearest neighbor is suppressed by $N^{-1/3}$, the nearest neighbor contributes an amount of order $G \bar{\rho} v$, which is comparable to the smooth-field jerk. This implies that the jerk remains stochastic and "grainy" even as $N \to \infty$ (at fixed total mass and density), whereas the acceleration becomes perfectly smooth in that limit because $a_{\mathrm{near}}/a_{\mathrm{smooth}} \to 0$.
