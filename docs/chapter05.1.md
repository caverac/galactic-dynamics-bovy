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
<!-- PROBLEM 5.5             -->
<!-- ======================= -->
## Problem 5.5

### Part a

In a razor-thin disk of radius $R$ with $N$ stars, the surface number density is $\Sigma_n = N/(\pi R^2)$. A star traversing the disk sweeps a strip of length $\sim 2R$ and width $2\,db$ (both sides) at impact parameter $b$:

$$
dN_{\mathrm{enc}} = \Sigma_n \times 2\,db \times 2R = \frac{4N}{\pi R}\,db.
$$

This is the key contrast with 3D, where the annular geometry gives $dN_{\mathrm{enc}} \propto b\,db$. The mean-squared velocity change per crossing is

$$
\begin{align}
\langle \Delta v^2 \rangle &= \int_{b_{\min}}^{R} dN_{\mathrm{enc}}\, (\delta v)^2 \\
&= \frac{16\,G^2 m^2 N}{\pi R\,v^2} \int_{b_{\min}}^{R} \frac{db}{b^2} \\
&= \frac{16\,G^2 m^2 N}{\pi R\,v^2} \left(\frac{1}{b_{\min}} - \frac{1}{R}\right) \\
&\approx \frac{16\,G^2 m^2 N}{\pi R\,v^2} \times \frac{1}{b_{\min}} \\
&= \frac{16\,G^2 m^2 N}{\pi R\,v^2} \times \frac{v^2}{2Gm} \\
&= \frac{8\,GmN}{\pi R}
\end{align}
$$

The number of crossings to relax is

$$
n_{\mathrm{cross}} = \frac{v^2}{\langle \Delta v^2 \rangle} = \frac{\pi\,R\,v^2}{8\,GmN}
$$

Applying the virial theorem $v^2 \sim GNm/R$:

$$
\frac{t_{\mathrm{relax}}}{t_{\mathrm{dyn}}} \sim \frac{\pi}{8} \sim \mathcal{O}(1)
$$

### Part b

The difference is purely geometric, in how encounters at impact parameter $b$ are counted:

- **3D (sphere):** encounters live in a cylindrical shell of area $2\pi b\,db$, so $dN_{\mathrm{enc}} \propto b\,db$. The factor of $b$ partially cancels the $1/b^2$ from $(\delta v)^2$, giving $\int db/b = \ln\Lambda$, all scales contribute roughly equally (logarithmically). Distant encounters collectively matter as much as close ones, and relaxation requires $\sim N/\ln N$ crossings.

- **2D (disk):** encounters live in a strip of width $2\,db$ with no extra factor of $b$, so $dN_{\mathrm{enc}} \propto db$. The integral becomes $\int db/b^2$, completely dominated by the closest encounters. Close encounters are so efficient that a single crossing already deflects a star by $\sim v$, making $t_{\mathrm{relax}} \sim t_{\mathrm{dyn}}$ regardless of $N$.

!!! warning "Incomplete"

    Parts c, d and e are missing.


<!-- ======================= -->
<!-- PROBLEM 5.6             -->
<!-- ======================= -->
## Problem 5.6


### Velocity dispersion
Start from the spherical Jeans equation for an isotropic system ($\beta = 0$):

$$
\frac{d(\nu\,\sigma^2)}{dr} = -\nu\,\frac{d\Phi}{dr}
$$

Integrate from $r$ to $\infty$, using the boundary condition $\nu\,\sigma^2 \to 0$ as $r \to \infty$:

$$
\nu(r)\,\sigma^2(r) = \int_r^\infty dr'\, \nu(r')\,\frac{d\Phi}{dr'}
$$

Now change the integration variable from $r'$ to $\Phi'$. In a spherical system $\Phi(r)$ is monotonic, so the substitution is one-to-one. The limits transform as $r' = r \Rightarrow \Phi' = \Phi(r)$ and $r' \to \infty \Rightarrow \Phi' \to 0$, giving

$$
\nu(r)\,\sigma^2(r) = \int_{\Phi(r)}^{0} d\Phi'\,\nu(\Phi')
$$

### Hernquist profile

The Hernquist profile has density and potential

$$
\rho(r) = \frac{M\,a}{2\pi\,r\,(r+a)^3}, \qquad
\Phi(r) = -\frac{GM}{r+a}
$$

Define $s = r/a$ and $q = -a\,\Phi/GM = 1/(1+s)$, inverting $s = (1-q)/q$, $1+s = 1/q$, so

$$
\rho(q) = \frac{M}{2\pi\,a^3}\,\frac{q^4}{1-q}
$$

Since $\Phi = -(GM/a)\,q$, we have $d\Phi = -(GM/a)\,dq$. The limits $\Phi(r) \to 0$ become $q_0 \to 0$ with $q_0 = 1/(1+s)$, so

$$
\begin{align}
\rho\,\sigma^2 &= \int_{\Phi(r)}^{0} \rho(\Phi')\,d\Phi'
= \frac{GM^2}{2\pi\,a^4}\int_0^{q_0}\frac{q^4}{1-q}\,dq \\
&= \frac{GM^2}{2\pi\,a^4}\int_0^{q_0}\left(-1 - q - q^2 - q^3 + \frac{1}{1-q}\right) dq
\end{align}
$$

So that

$$
\sigma^2(s) = \frac{GM}{a}\,s(1+s)^3\left[
\ln\!\left(\frac{1+s}{s}\right)
- \frac{1}{1+s}
- \frac{1}{2(1+s)^2}
- \frac{1}{3(1+s)^3}
- \frac{1}{4(1+s)^4}
\right]
$$

<!-- ======================= -->
<!-- PROBLEM 5.7             -->
<!-- ======================= -->
## Problem 5.7 🌶️

Write $f(\mathcal{E}, L) = L^{-2\beta}\,g(\mathcal{E})$ where

$$
g(\mathcal{E}) \propto \tilde{\mathcal{E}}^{\,5/2-\beta}\; {}_2F_1\!\bigl(5-2\beta,\; 1-2\beta;\; 7/2-\beta;\; \tilde{\mathcal{E}}\bigr).
$$

The velocity distribution at radius $r$ is

$$
\begin{align}
p(v \mid r) &\propto v^2 \int d\Omega_v\, f(\mathcal{E}, L) \\
&= v^2 g(\mathcal{E}) (rv)^{-2\beta} \int d\Omega_v\, (\sin\theta)^{-2\beta} \\
&\propto v^{2-2\beta}\; g(\mathcal{E})
\end{align}
$$

**Concentration on circular orbits.** The tangential-orbit limit corresponds to $\beta \to -\infty$. In this regime, $f(\mathcal{E}, L) \propto L^{2|\beta|}\, g(\mathcal{E})$. Introduce the circularity

$$
\eta \equiv \frac{L}{L_c(\mathcal{E})},
$$

where $L_c(\mathcal{E})$ is the angular momentum of the circular orbit at energy $\mathcal{E}$. By definition $0 \le \eta \le 1$,  with equality $\eta = 1$ only for circular orbits. Thus,

$$
f(\mathcal{E}, L) \propto \eta^{2|\beta|}\, L_c(\mathcal{E})^{2|\beta|}\, g(\mathcal{E}).
$$

At fixed radius $r$ and speed $v$, the angular momentum is

$$
L = r v \sin\theta,
$$

so

$$
\eta = \frac{r v \sin\theta}{L_c(\mathcal{E})} \le \frac{r v}{L_c(\mathcal{E})} \equiv Q_r(v) \le 1.
$$

The upper bound $Q_r(v)=1$ is reached only when:

- the velocity is purely tangential ($\sin\theta = 1$),
- and $v = v_c(r)$ is the circular speed at radius $r$.

After integrating over velocity directions, the speed distribution behaves as

$$
p(v \mid r) \propto \left[Q_r(v)\right]^{2|\beta|}
$$

Since:

- $Q_r(v) < 1$ for all $v \ne v_c(r)$,
- $Q_r(v) = 1$ only at $v = v_c(r)$,

it follows that

$$
\left[Q_r(v)\right]^{2|\beta|} \to 0 \quad \text{for } v \ne v_c(r),
$$

as $|\beta| \to \infty$.

Therefore, the speed distribution collapses to a delta function:

$$
p(v \mid r) \xrightarrow[\beta \to -\infty]{} \delta\!\bigl(v - v_c(r)\bigr).
$$

<!-- ======================= -->
<!-- REFERENCES              -->
<!-- ======================= -->
## References
\bibliography
