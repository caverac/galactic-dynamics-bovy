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

<!-- ======================= -->
<!-- PROBLEM 5.4             -->
<!-- ======================= -->
## Problem 5.4

### Part a

The mean-squared velocity change per crossing comes from two-body encounters with both stars and dark-matter particles. For a species with individual mass $m_i$ and total number $N_i$:

$$
\langle \Delta v^2 \rangle_i = 8 N_i \left(\frac{G m_i}{R v}\right)^2 \ln \Lambda_i
$$

Writing $N_i = M_i / m_i$ where $M_i$ is the total mass of species $i$:

$$
\langle \Delta v^2 \rangle_i \propto M_i\, m_i
$$

The ratio of dark-matter to stellar contributions is therefore

$$
\frac{\langle \Delta v^2 \rangle_{\mathrm{DM}}}{\langle \Delta v^2 \rangle_*} = \frac{M_{\mathrm{DM}}\, m_{\mathrm{DM}}}{M_*\, m} = \frac{f}{1-f}\,\frac{m_{\mathrm{DM}}}{m} \tag{5.4.1}
$$

The mass-fraction ratio $f/(1-f)$ is of order unity (e.g. $\approx 1$ for $f = 0.5$), but $m_{\mathrm{DM}}/m \ll 1$ so the product is negligible. That is, if dark matter is made of particles, then its discreteness negligibly contributes to two-body relaxation compared to the stars.

With total galaxy mass $M$, dark-matter fraction $f$, stellar mass $M_* = (1-f)M$, and $N_* = (1-f)M/m$ stars, the dominant contribution is

$$
\langle \Delta v^2 \rangle \approx 8 N_* \left(\frac{Gm}{Rv}\right)^2 \ln \Lambda
$$

The velocity dispersion is still set by the total mass via the virial theorem, $v^2 \sim GM/R$, so the relaxation time is

$$
t_{\mathrm{relax}} = \frac{v^2}{\langle \Delta v^2 \rangle / t_{\mathrm{cross}}} = \frac{M}{8(1-f)\,m\,\ln\Lambda}\,t_{\mathrm{cross}}
$$

With the Coulomb logarithm $\Lambda \approx R v^2/(Gm) = M/m$, this becomes

$$
t_{\mathrm{relax}} = \frac{1}{1-f}\,\frac{M/m}{8\ln(M/m)}\,t_{\mathrm{cross}}
$$

Compared to a galaxy of the same total mass with no dark matter, the relaxation time is longer by a factor of $1/(1-f)$. Dark matter contributes to the smooth gravitational potential (setting $v$) but, because $m_{\mathrm{DM}} \ll m$, it does not contribute to the graininess that drives two-body relaxation. Fewer stars for the same total mass means a smoother system and slower relaxation.

### Part b

From Eq. (5.4.1), dark-matter encounters become important when the ratio is of order unity:

$$
\frac{f}{1-f}\,\frac{m_{\mathrm{DM}}}{m} \sim 1
\quad \Longrightarrow \quad
m_{\mathrm{DM}} \sim \frac{1-f}{f}\,m
$$

For $f \approx 0.5$ this gives $m_{\mathrm{DM}} \sim m \sim M_\odot$. This is precisely the regime of MACHOs [massive compact halo objects; @griest1991] or primordial black holes of stellar mass. For any $m_{\mathrm{DM}} \gtrsim \frac{1-f}{f}\,m$, the dark-matter particles contribute significantly to two-body relaxation and must be included in the calculation.


### Part c

With $m_{\mathrm{DM}} = 100\,M_\odot$ and $f = 0.5$, Eq. (5.4.1) gives a ratio of $100$: dark-matter encounters now dominate relaxation. Including both species, the general relaxation time is

$$
\begin{align}
t_{\mathrm{relax}} &= \frac{M^2}{8\bigl[(1-f)\,M\,m\,\ln\Lambda_* + f\,M\,m_{\mathrm{DM}}\,\ln\Lambda_{\mathrm{DM}}\bigr]}\,t_{\mathrm{cross}} \\
&= \frac{M}{8\bigl[(1-f)\,m\,\ln\Lambda_* + f\,m_{\mathrm{DM}}\,\ln\Lambda_{\mathrm{DM}}\bigr]}\,t_{\mathrm{cross}}
\end{align}
$$

where $\Lambda_i \approx M/m_i$. For the Milky Way within the visible region we adopt $M \approx 10^{11}\,M_\odot$, $R \approx 10\,\mathrm{kpc}$, $v \approx 220\,\mathrm{km/s}$, with $m = 1\,M_\odot$:

$$
\Lambda_* = M/m = 10^{11} \;\Rightarrow\; \ln\Lambda_* \approx 25.3
\qquad
\Lambda_{\mathrm{DM}} = M/m_{\mathrm{DM}} = 10^{9} \;\Rightarrow\; \ln\Lambda_{\mathrm{DM}} \approx 20.7
$$

The crossing time is

$$
t_{\mathrm{cross}} = \frac{2\pi R}{v} = \frac{2\pi \times 10\,\mathrm{kpc}}{220\,\mathrm{km/s}} \approx 0.28\,\mathrm{Gyr}
$$

And

$$
(1-f)\,m\,\ln\Lambda_* + f\,m_{\mathrm{DM}}\,\ln\Lambda_{\mathrm{DM}} \approx 1048\,M_\odot
$$

confirming that the primordial black holes contribute $\sim\!99\%$ of the relaxation. Thus

$$
t_{\mathrm{relax}}
\approx \frac{10^{11}}{8 \times 1048} \times 0.28\,\mathrm{Gyr}
\approx 3.3 \times 10^{6}\,\mathrm{Gyr}
$$

For comparison, a pure-stellar Milky Way ($f = 0$, $N = 10^{11}$) would have $t_{\mathrm{relax}} \approx 1.4 \times 10^{8}\,\mathrm{Gyr}$, about 40 times longer. Even so, $3.3 \times 10^{6}\,\mathrm{Gyr} \gg t_{\mathrm{Hubble}} \approx 14\,\mathrm{Gyr}$, so the Milky Way remains safely collisionless even in this extreme scenario.

### Part d

The dark-matter mass enclosed within the cluster's half-light radius is

$$
M_{\mathrm{DM}} = \rho_{\mathrm{DM}}\,\tfrac{4}{3}\pi\,r_h^3
= 0.3\,\frac{M_\odot}{\mathrm{pc}^3} \times \tfrac{4}{3}\pi\,(13\,\mathrm{pc})^3
\approx 2{,}760\,M_\odot
$$

giving $N_{\mathrm{DM}} = M_{\mathrm{DM}}/m_{\mathrm{DM}} \approx 28$ primordial black holes inside the cluster, alongside $N_* = 5{,}000$ stars. The Coulomb logarithms for each species are

$$
\Lambda_i = \frac{r_h\,\sigma^2}{G\,m_i}
\quad\Rightarrow\quad
\ln\Lambda_* \approx 11.2,\quad \ln\Lambda_{\mathrm{DM}} \approx 6.6
$$

The contributions to $\langle \Delta v^2 \rangle$ per crossing are proportional to $N_i\,m_i^2\,\ln\Lambda_i$:

$$
N_*\,m^2\,\ln\Lambda_* = 5{,}000 \times 1 \times 11.2 = 5.6 \times 10^4\,M_\odot^2
$$

$$
N_{\mathrm{DM}}\,m_{\mathrm{DM}}^2\,\ln\Lambda_{\mathrm{DM}} = 28 \times 10^4 \times 6.6 = 1.85 \times 10^6\,M_\odot^2
$$

so the 28 primordial black holes contribute $\sim\!97\%$ of the relaxation despite being vastly outnumbered. The number of crossings to relax is

$$
n_{\mathrm{cross}}
= \frac{\sigma^4\,r_h^2}{8\,G^2\,\sum_i N_i\,m_i^2\,\ln\Lambda_i}
\approx 375
$$

With $t_{\mathrm{cross}} = 2\pi\,r_h/\sigma \approx 16\,\mathrm{Myr}$:

$$
t_{\mathrm{relax}} \approx 375 \times 16\,\mathrm{Myr} \approx 6\,\mathrm{Gyr}
$$

This is less than the age of the universe, so $100\,M_\odot$ primordial black holes could plausibly destroy the Eridanus II star cluster through two-body encounters. For comparison, without dark matter the cluster would have $t_{\mathrm{relax}} \approx 200\,\mathrm{Gyr}$, safely collisionless.

<!-- ======================= -->
<!-- REFERENCES              -->
<!-- ======================= -->
## References
\bibliography
