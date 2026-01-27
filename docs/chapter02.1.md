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