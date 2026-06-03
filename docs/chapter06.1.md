# Chapter 6: Masses of spherical stellar systems

<!-- ======================= -->
<!-- PROBLEM 6.1             -->
<!-- ======================= -->
## Problem 6.1

From Eqns. (6.12) and (2.24)

$$
M(<r) = -\frac{r\,\sigma_r^2}{G}\left(\frac{d\ln\nu}{d\ln r} + \frac{d\ln\sigma_r^2}{d\ln r} + 2\beta\right).
$$

For the Milky Way halo BHB sample [@xue2008], the tracer number density falls as $\nu \propto r^{-\gamma}$ with $\gamma = 3.5$, so $d\ln\nu/d\ln r = -3.5$. Taking $\sigma_r \approx$ const removes the dispersion term, $d\ln\sigma_r^2/d\ln r = 0$, leaving

$$
M(<r) = \frac{r\,\sigma_r^2}{G}\,(\gamma - 2\beta) = \frac{r\,\sigma_r^2}{G}\,(3.5 - 2\beta).
$$

And from here

$$
\frac{dM}{dr} = \frac{\sigma_r^2}{G}\,(3.5 - 2\beta).
$$

For the enclosed mass not to increase we need $dM/dr \le 0$, i.e.

$$
\gamma - 2\beta \le 0 \quad\Longrightarrow\quad \beta \ge \frac{\gamma}{2} = 1.75.
$$

Which is not possible since $\beta$ is bounded above by $\beta = 1$.

<!-- ======================= -->
<!-- PROBLEM 6.2             -->
<!-- ======================= -->
## Problem 6.2

$$
M_{1/2} = \frac{4\,\sigma_{\mathrm{los}}^2\,R_e}{G}
\approx 930\,\left(\frac{\sigma_{\mathrm{los}}}{\mathrm{km\,s^{-1}}}\right)^2\left(\frac{R_e}{\mathrm{pc}}\right)\,M_\odot,
$$

using $G = 4.30\times10^{-3}\,\mathrm{pc}\,(\mathrm{km/s})^2/M_\odot$. With $\sigma_{\mathrm{los}} \approx 4\,\mathrm{km\,s^{-1}}$ and $R_e \approx 5\,\mathrm{pc}$,

$$
M_{1/2} \approx 930 \times 4^2 \times 5 \approx 7.4\times10^{4}\,M_\odot,
$$

and the total mass is about twice this (the half-light radius encloses half the tracers),

$$
M_{\mathrm{tot}} \approx 2\,M_{1/2} \approx 1.5\times10^{5}\,M_\odot.
$$

Milky Way globular clusters have half-light radii of a few pc and masses $\sim10^4$-$10^6\,M_\odot$ [@harris1996], which puts the $R_e \approx 5\,\mathrm{pc}$ in the GC regime. Dwarf galaxies are far more extended — Local Group dwarfs have $R_e \gtrsim 100\,\mathrm{pc}$, with even the smallest ultra-faints rarely below $\sim30\,\mathrm{pc}$ [@mcconnachie2012], and dwarf spheroidals appear to have a characteristic size floor near $\sim100$-$120\,\mathrm{pc}$ [@gilmore2007]. A few-pc system therefore cannot be a dwarf galaxy.
