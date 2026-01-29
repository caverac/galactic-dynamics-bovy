<!-- ======================= -->
<!-- PROBLEM 2.11             -->
<!-- ======================= -->
## Problem 2.11

Consider a potential of the form

$$
\Phi(r) = -\frac{GM_p}{a}\frac{1}{(1 + r^p/a^p)^{1/p}}
$$


### Density
The density can be directly calculated from Poisson's equation.

$$
\begin{align}
\rho(r) &= \frac{1}{4\pi G}\nabla^2\Phi(r) \\
&= \frac{1}{4\pi G}\frac{1}{r^2}\frac{d}{dr}\left(r^2\frac{d\Phi}{dr}\right) \\
&= \frac{M_p}{4\pi}\frac{1}{r^2}\frac{d}{dr}\left(\frac{r^2 (r/a)^{p-1}}{a^2(1 + (r/a)^p)^{1 + 1/p}}\right) \\
&= \frac{M_p(p+1)}{4\pi a^3}\frac{(r/a)^{p-2}}{(1 + (r/a)^p)^{2 + 1/p}} \\
\end{align}
$$

### Mass

The mass can also be directly calculated from the potential

$$
\begin{align}
M(r) &= \frac{1}{G}r^2\frac{d\Phi}{dr} \\
&= \frac{M_p r^2 (r/a)^{p-1}}{a^2(1 + (r/a)^p)^{1 + 1/p}} \\
&= M_p\frac{(r/a)^{p+1}}{(1 + (r/a)^p)^{1 + 1/p}} \\
\end{align}
$$

Total mass is $M_p$ as $r\to\infty$, so the $r_{1/2}$ can be found by solving

$$
M(r_{1/2}) = \frac{M_p}{2} = M_p\frac{(r_{1/2}/a)^{p+1}}{(1 + (r_{1/2}/a)^p)^{1 + 1/p}}
$$

Define $x = (r_{1/2}/a)^p$, so that the equation becomes

$$
\frac{1}{2} = \frac{x^{(p+1)/p}}{(1 + x)^{1 + 1/p}} = \left(\frac{x}{1 + x}\right)^{(p+1)/p}
$$

Threfore

$$
r_{1/2} = a\left(2^{p/(p+1)} - 1\right)^{-1/p}
$$

<!-- ======================= -->
<!-- PROBLEM 2.12             -->
<!-- ======================= -->
## Problem 2.12

Instead of calculating the projected mass using $\Sigma(R) = 2\int_0^\infty dz\, \rho(\sqrt{R^2 + z^2})$ I will use the Abel form

$$
\Sigma(R) = 2\int_R^\infty dr\,\frac{r\rho(r)}{\sqrt{r^2 - R^2}},
$$

with

$$
\rho(r) = \frac{(3 - \gamma)M}{4\pi a^3} \frac{1}{y^\gamma(1 + y)^{4 - \gamma}},
$$

Defining $x = R/a$ and $y = r/a$, this becomes

$$
\Sigma(R) = \frac{M}{a^2}s(x), \quad s(x) = \frac{(3 - \gamma)}{2\pi}\int_x^\infty dy\,\frac{y^{1 - \gamma}}{(1 + y)^{4 - \gamma}\sqrt{y^2 - x^2}}. \tag{2.12.1}
$$

and the projected mass

$$
M_{\mathrm{proj}}(<R) = M m_p(x), \quad m_p(x) = 2\pi\int_0^x dx'\, x' s(x'). \tag{2.12.2}
$$

The effective radius requires to find $x_e$ such that $m_p(x_e) = 1/2$. Also, the half-mass radius in 3D is given by $y_{1/2} = (2^{1/(3 - \gamma)} - 1)^{-1}$.

Before I go into the details of the numerical implementation, I will plot the effective radius as a function of $\gamma$.

```python
from galactic_dynamics_bovy.chapter02.effective_radius_gamma import plot_effective_radius_gamma
plot_effective_radius_gamma()
```

![Effective radius Gamma model](assets/generated/effective_radius_gamma.png)

*Figure 2.12: Effective radius as a function of $\gamma$. For large $\gamma$ the mass becomes more centrally concentrated, which makes numerically convergence a bit more challenging. Line $R_e/r_{1/2} = 3/4$ added for reference.*

Now, the numerical considerations

1. I pre-compute $s(x)$ on a grid in $x$ for a given $\gamma$ using `scipy.integrate.quad` to evaluate the integral in equation (2.12.1). Same as $m_p(x)$.
2. The values are interpolated using `scipy.interpolate.interp1d`. If a value outside the grid is requested, I directly use equations (2.12.1) and (2.12.2) to compute the values.
3. Instead of optimizing for $x_e$ such that $m_p(x_e) = 1/2$, I just evalute the interpolation $x = x(m_p)$
4. All these interpolations and numerical integrations deal with numerical divergences at small $x$, unfortunately large values of $\gamma$ make this more challenging, as the density becomes more centrally concentrated. I use logarithmic grids and careful handling of small values to mitigate this, however numerical noise remains for large $\gamma$.
