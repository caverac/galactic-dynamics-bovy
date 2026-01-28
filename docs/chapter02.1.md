# Chapter 2: Gravitation

<!-- ======================= -->
<!-- PROBLEM 2.1             -->
<!-- ======================= -->
## Problem 2.1

For a homogeneous sphere $v_c \sim r$

$$
v_c(r) = \sqrt{\frac{4\pi G \rho_0}{3}} r
$$

![IC 2574 Rotation Curve](assets/generated/ic2574_rotation_curve.png)

*Figure 2.1: Data from the THINGS survey [@deblok2008] shows that the circular velocity rises approximately linearly with radius for $r < 6$ kpc. Linear fit shown in dashed line.*

$$v_c \approx 6.75 \cdot r + 8.18 \text{ km/s}$$

with $R^2 = 0.949$, from this we can estimate the central density of IC 2574 as

$$
\rho_0 = \frac{3}{4 \pi G} \left(\frac{v_c}{r}\right)^2 \approx 1.27 \times 10^7 M_\odot \text{ kpc}^{-3}
$$

Script to generate figure:

```python
from galactic_dynamics_bovy.chapter02.ic_2574 import plot_rotation_curve
plot_rotation_curve()
```

<!-- ======================= -->
<!-- PROBLEM 2.2             -->
<!-- ======================= -->
## Problem 2.2

### Numerical solution

The connection between virial radius $r_{\mathrm{vir}}$, virial mass $M_{\mathrm{vir}}$, and overdensity $\Delta_v$ is

$$
\frac{3\rho_0}{c^3}f(c) = \Delta_v \rho_{\mathrm{crit}} \quad f(c) = \ln(1+c) - \frac{c}{1+c} \quad r_{\mathrm{vir}} = c a,
$$

or equivalently

$$
\frac{3f(c)}{c^3} = y \quad y = \frac{\Delta_v \rho_{\mathrm{crit}}}{\rho_0}.
$$

The expression $r_{\mathrm{vir}} = a c(y)$ then gives the virial radius as a function of overdensity $\Delta_v$. However, there is no closed-form solution for $c(y)$, so we must solve numerically.

![NFW Virial Mass and Radius](assets/generated/nfw_rvir_mvir_delta.png)

*Figure 2.2: Virial mass (top) and virial radius (bottom) as a function of overdensity $\Delta_v$ for an NFW halo with $\rho_0 = 3.5 \times 10^6\,M_\odot/\mathrm{kpc}^3$ and scale radius $a = 16$ kpc. At $\Delta_v = 200$, this halo has $r_{\mathrm{vir}} \approx 95$ kpc and $M_{\mathrm{vir}} \approx 2 \times 10^{11}\,M_\odot$.*

Script to generate figure:

```python
from galactic_dynamics_bovy.chapter02.rvir_mvir import plot_rvir_mvir_delta
plot_rvir_mvir_delta()
```

That being said, there are reasonable approximations for $c(y)$ at high and low $y$:

### Small overdensity

In this regime

$$
\begin{align}
f(c) &= \ln(1+c) - \frac{c}{1+c} \\
&= \left(c - \frac{c^2}{2} + \frac{c^3}{3} - \ldots\right) - \left(c - c^2 + c^3 - \ldots\right) \\
&= \frac{c^2}{2} - \frac{2 c^3}{3} + \mathcal{O}(c^4)
\end{align}
$$

such that

$$
\frac{2f(c)}{c^3} = \frac{3}{2c} - 2 + \mathcal{O}(c).
$$

To leading order

$$
c \approx \frac{3}{2y},
$$

or equivalently

$$
r_{\mathrm{vir}} \approx \frac{3a}{2y} = \frac{3a \rho_0}{2 \Delta_v \rho_{\mathrm{crit}}} \sim \Delta_v^{-1}
$$

### Large overdensity
In this regime

$$
\begin{align}
f(c) &= \ln(1+c) - \frac{c}{1+c} \\
&= \ln c - 1 + \mathcal{O}(1/c)
\end{align}
$$

therefore

$$
c^3 \sim \frac{2\ln c}{y} \implies c \sim \left(\frac{2 \ln c}{y}\right)^{1/3}
$$

Up to a logarithmic correction, we have

$$
r_{\mathrm{vir}} \sim c \sim y^{-1/3} \sim \Delta_v^{-1/3}
$$

In summary, we have

$$
r_{\mathrm{vir}} \sim \begin{cases}\Delta_v^{-1}, & \Delta_v \ll 1 \\ \Delta_v^{-1/3}, & \Delta_v \gg 1\end{cases}
$$

### Intuition behind the scaling

If $\rho(r) \sim r^{-\gamma}$ then $M(<r) \sim r^{3-\gamma}$, and $\overline{\rho}(<r) \sim r^{-\gamma}$. Setting $\overline{\rho}(<r) = \Delta_v \rho_{\mathrm{crit}}$ gives the scaling of $r_{\mathrm{vir}}$ with $\Delta_v$.

$$
r_{\mathrm{vir}} \sim \Delta_v^{-1/\gamma}
$$

The NFW profile has $\gamma \approx 1$ at small $r$ and $\gamma \approx 3$ at large $r$, giving the scalings derived above. The reason for the scaling is then the two-power nature of the NFW profile.

!!! note "Milky Way NFW Parameters"

    Recent estimates of the Milky Way's virial mass span a wide range. Using Gaia DR3 data, [@wang2023] found $M_{\mathrm{vir}} = (6.5 \pm 0.3) \times 10^{11}\,M_\odot$ with concentration $c \approx 14.5$ for an NFW model. Earlier work by [@eilers2019] estimated $M_{\mathrm{vir}} = (7.25 \pm 0.25) \times 10^{11}\,M_\odot$. Studies using globular clusters [@posti2019] suggest a higher mass of $M_{\mathrm{vir}} = (1.3 \pm 0.3) \times 10^{12}\,M_\odot$.


<!-- ======================= -->
<!-- REFERENCES              -->
<!-- ======================= -->
## References
\bibliography
