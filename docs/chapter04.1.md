# Chapter 4: Orbits in spherical mass distributions

<!-- ======================= -->
<!-- PROBLEM 4.1             -->
<!-- ======================= -->
## Problem 4.1

In a homogeneous spherical potential, the radial frequency is constant, independent of energy and angular momentum, $T_r = \pi/\omega$.

<!-- ======================= -->
<!-- PROBLEM 4.2             -->
<!-- ======================= -->
## Problem 4.2

### Qualitative argument


There are two limiting cases for the ratio of the radial period $T_r$ to the azimuthal period $T_\psi$:

| Potential | Rotation curve | $T_r / T_\psi$ |
|---|---|---|
| Homogeneous sphere | $v_c \propto r$ | $1/2$ |
| Kepler | $v_c \propto r^{-1/2}$ | $1$ |

with the isochrone interpolating between them: $1/2 \leq T_r/T_\psi \leq 1$.

The Milky Way has a roughly flat rotation curve near the Sun ($v_c \approx \mathrm{const}$), which lies between the rising rotation curve of the homogeneous sphere and the falling Keplerian one. Therefore $T_r/T_\psi$ should be somewhere between $1/2$ and $1$, giving

$$
T_r \sim 110\text{-}220\,\mathrm{Myr}.
$$

### More refined argument

For a power-law rotation curve $v_c \propto r^\beta$, the circular frequency is $\Omega = v_c/r \propto r^{\beta - 1}$, so $\Omega^2 \propto r^{2(\beta-1)}$. The epicyclic frequency is

$$
\kappa^2 = R\frac{d\Omega^2}{dR} + 4\Omega^2 = 2(\beta - 1)\Omega^2 + 4\Omega^2 = 2(\beta + 1)\,\Omega^2,
$$

giving

$$
\frac{T_r}{T_\psi} = \frac{\Omega}{\kappa} = \frac{1}{\sqrt{2(\beta + 1)}}.
$$

This reproduces the two limits: $\beta = 1 \to T_r/T_\psi = 1/2$ and $\beta = -1/2 \to T_r/T_\psi = 1$.

For a flat rotation curve ($\beta = 0$):

$$
\frac{T_r}{T_\psi} = \frac{1}{\sqrt{2}},
$$

and with $T_\psi \approx 220\,\mathrm{Myr}$:

$$
T_r = \frac{220\,\mathrm{Myr}}{\sqrt{2}} \approx 155\,\mathrm{Myr}.
$$

<!-- ======================= -->
<!-- PROBLEM 4.3             -->
<!-- ======================= -->
## Problem 4.3

### Part a

With

$$
s = 1 + \sqrt{1 + \frac{r^2}{b^2}},
$$

the isochrone potential can be written as

$$
\begin{align}
\Phi(r) &= -\frac{GM}{b + \sqrt{b^2 + r^2}} \\
&= -\frac{GM}{b \left(1 + \sqrt{1 + r^2/b^2}\right)} \\
&= -\frac{GM}{b s}.
\end{align}
$$

We can invert $r$ to get

$$
r^2 = b^2(s^2 - 2s) = b^2 s(s-2) = f(s).
$$

### Part b

Becasue we are going to need it later, call

$$
g(r) = 2[E - \Phi(r)] - \frac{L^2}{r^2} = 2E +2\frac{GM}{b s} - \frac{L^2}{b^2 s(s-2)},
$$

So that

$$
s(s - 2) g(r) = 2Es^2 - 2\left(2E - \frac{GM}{b}\right)s - \left(\frac{4GM}{b} + \frac{L^2}{b^2}\right). \tag{4.3.1}
$$

The turning points are solutions to the equation $g(r) = 0$,

$$
2Es^2 - 2\left(2E - \frac{GM}{b}\right)s - \left(\frac{4GM}{b} + \frac{L^2}{b^2}\right) = 0,
$$

whose solutions are

$$
4Es_{1,2} = 2\left(2E - \frac{GM}{b}\right) \pm \sqrt{4\left(2E - \frac{GM}{b}\right)^2 + 8E\left(\frac{4GM}{b} + \frac{L^2}{b^2}\right)}.
$$

Therefore

$$
s(s - 2) g(r) = 2E(s - s_1)(s - s_2) = -2E(s_2 - s)(s - s_1).
$$

### Part c

$$
\begin{align}
T_r &= 2\int_{r_p}^{r_a} dr\, \frac{1}{\sqrt{g(r)}} \\
&= 2\int_{s_1}^{s_2} ds\,\frac{dr}{ds} \frac{1}{\sqrt{g(r)}} \\
&= 2\int_{s_1}^{s_2} ds\, \frac{b(s - 1)}{\sqrt{s(s-2)g(r)}} \\
&= \frac{2b}{\sqrt{-2E}} \int_{s_1}^{s_2} ds\, \frac{s - 1}{\sqrt{(s_2 - s)(s - s_1)}} \\
&= \frac{2b}{\sqrt{-2E}} \frac{\pi}{2}(s_1 + s_2 - 2) \\
&\stackrel{(4.3.1)}{=} \frac{\pi b}{\sqrt{-2E}}\left[\frac{4}{4E}\left(2E - \frac{GM}{b} \right) - 2\right] \\
&= -\frac{\pi b}{\sqrt{-2E}}\frac{GM}{bE} \\
&= \frac{2\pi GM}{(-2E)^{3/2}}.
\end{align}
$$

### Part d


$$
\begin{align}
\Delta\psi &= 2\int_{r_p}^{r_a} \frac{L}{r^2} \frac{dr}{\sqrt{g(r)}} \\
&= 2\int_{s_1}^{s_2} ds\, \frac{L}{b^2 s(s-2)} \cdot \frac{b(s-1)}{\sqrt{s(s-2)}} \cdot \frac{\sqrt{s(s-2)}}{\sqrt{-2E}\sqrt{(s_2-s)(s-s_1)}} \\
&= \frac{2L}{b\sqrt{-2E}} \int_{s_1}^{s_2} ds\, \frac{s - 1}{s(s-2)\sqrt{(s_2 - s)(s - s_1)}} \\
&= \Delta\psi = \frac{\pi L}{b\sqrt{-2E}}\left(\frac{1}{\sqrt{s_1 s_2}} + \frac{1}{\sqrt{(s_1-2)(s_2-2)}}\right).
\end{align}
$$



From Vieta's formulas applied to the quadratic (4.3.1),

$$
s_1 s_2 = -\frac{1}{2E}\left(\frac{4GM}{b} + \frac{L^2}{b^2}\right) = \frac{4GMb + L^2}{2b^2(-E)},
$$

and

$$
(s_1-2)(s_2-2) = s_1 s_2 - 2(s_1 + s_2) + 4 = \frac{L^2}{2b^2(-E)},
$$

where in the last step we used $s_1 + s_2 = 2 - GM/(bE)$. Therefore

$$
\frac{1}{\sqrt{s_1 s_2}} = \frac{b\sqrt{-2E}}{\sqrt{L^2 + 4GMb}}, \qquad \frac{1}{\sqrt{(s_1-2)(s_2-2)}} = \frac{b\sqrt{-2E}}{|L|}.
$$

Substituting back,

$$
\Delta\psi = \frac{\pi L}{b\sqrt{-2E}}\left(\frac{b\sqrt{-2E}}{\sqrt{L^2 + 4GMb}} + \frac{b\sqrt{-2E}}{|L|}\right) = \pi\,\mathrm{sign}(L)\left(1 + \frac{|L|}{\sqrt{L^2 + 4GMb}}\right).
$$

<!-- ======================= -->
<!-- PROBLEM 4.4             -->
<!-- ======================= -->
## Problem 4.4

The latitudinal action is

$$
J_\theta = \frac{1}{2\pi}\oint d\theta\, p_\theta
$$

From [Problem 3.3](chapter03.1.md#problem-33), the conjugate momenta per unit mass are $p_r = \dot{r}$, $p_\theta = r^2\dot{\theta}$, and $p_\phi = r^2\sin^2\theta\,\dot{\phi} = L_z$. Since only the tangential components of $\mathbf{v}$ contribute to $\mathbf{L} = \mathbf{r}\times\mathbf{v}$:

$$
L^2 = r^4\dot{\theta}^2 + r^4\sin^2\theta\,\dot{\phi}^2 = p_\theta^2 + \frac{L_z^2}{\sin^2\theta},
$$

so the conjugate momentum is

$$
p_\theta = \sqrt{L^2 - \frac{L_z^2}{\sin^2\theta}}.
$$

The orbit in $\theta$ oscillates between turning points $\theta_{\min}$ and $\theta_{\max} = \pi - \theta_{\min}$, where $p_\theta = 0$ gives $\sin\theta_{\min} = |L_z|/L$. Therefore

$$
J_\theta = \frac{1}{\pi}\int_{\theta_{\min}}^{\pi - \theta_{\min}} d\theta\,\sqrt{L^2 - \frac{L_z^2}{\sin^2\theta}}.
$$

Substituting $u = \cos\theta$ (so $du = -\sin\theta\,d\theta$) with $u_0 = \cos\theta_{\min} = \sqrt{1 - L_z^2/L^2}$:

$$
J_\theta = \frac{L}{\pi}\int_{-u_0}^{u_0} du\,\frac{\sqrt{u_0^2 - u^2}}{1 - u^2},
$$

where I used $L^2 - L^2 u^2 - L_z^2 = L^2(u_0^2 - u^2)$.

Now substitute $u = u_0\sin t$ (so $du = u_0\cos t\,dt$):

$$
J_\theta = \frac{2Lu_0^2}{\pi}\int_0^{\pi/2} dt\,\frac{\cos^2 t}{1 - u_0^2\sin^2 t}.
$$

Writing $\cos^2 t = 1 - \sin^2 t$ and using

$$
\frac{\sin^2 t}{1 - u_0^2\sin^2 t} = \frac{1}{u_0^2}\left(\frac{1}{1 - u_0^2\sin^2 t} - 1\right),
$$

the integrand simplifies to

$$
\frac{\cos^2 t}{1 - u_0^2\sin^2 t} = \frac{u_0^2 - 1}{u_0^2}\cdot\frac{1}{1 - u_0^2\sin^2 t} + \frac{1}{u_0^2}.
$$

Using the standard result $\displaystyle\int_0^{\pi/2}\frac{dt}{1 - a^2\sin^2 t} = \frac{\pi}{2\sqrt{1-a^2}}$ for $a^2 < 1$, and noting $\sqrt{1 - u_0^2} = |L_z|/L$:

$$
\int_0^{\pi/2}\frac{dt}{1 - u_0^2\sin^2 t} = \frac{\pi}{2|L_z|/L} = \frac{\pi L}{2|L_z|}.
$$

Putting it all together:

$$
\begin{align}
J_\theta &= \frac{2L}{\pi}\left[\left(u_0^2 - 1\right)\cdot\frac{\pi L}{2|L_z|} + \frac{\pi}{2}\right] \\
&= \frac{2L}{\pi}\left[-\frac{L_z^2}{L^2}\cdot\frac{\pi L}{2|L_z|} + \frac{\pi}{2}\right] \\
&= \frac{2L}{\pi}\cdot\frac{\pi}{2}\left[1 - \frac{|L_z|}{L}\right] \\
&= L - |L_z|.
\end{align}
$$

<!-- ======================= -->
<!-- PROBLEM 4.5             -->
<!-- ======================= -->
## Problem 4.5

For large $b$

$$
\sqrt{b^2 + r^2} = b\sqrt{1 + \frac{r^2}{b^2}} = b\left(1 + \frac{1}{2}\frac{r^2}{b^2} - O(b^{-4}\right),
$$

so

$$
b + \sqrt{b^2 + r^2} = 2b + \frac{r^2}{2b} - O(b^{-3}),
$$

and the potential becomes

$$
\begin{align}
\Phi(r) &= -\frac{GM}{b + \sqrt{b^2 + r^2}} \\
&= -\frac{GM}{2b}\left(1 - \frac{r^2}{4b^2} + O(b^{-4})\right) \\
&= -\frac{GM}{2b} + \frac{GM}{8b^3}r^2 + O(b^{-5}).
\end{align}
$$

Up to an additive constant, this is a homogeneous sphere potential with $\omega^2 = GM/(4b^3)$. But the key point is that we must keep $GM/b^3$ fixed as $b \to \infty$ to get a non-trivial limit, in other words $M$ must scale as $b^3$ for the isochrone to approach a homogeneous sphere.

The "paradox" arises if we try to keep $M$ fixed as $b \to \infty$, the potential becomes trivial $\Phi(r) \to 0$, bound orbits disapper and trying to reason about perior becomes a meaningless exercise.

Note that $\Phi(0) = -GM/(2b) \to 0$ as $b \to \infty$, the depth of the potential well goes to zero if $M$ is fixed, abd tge bound energy windoe squeezes to zero, so holding $E$ fixed as $b \to \infty$ also leads to a trivial limit. The comparison is to hold the oscillator energy above the minimum $E - \Phi(0) = E + GM/(2b)$ fixed as $b \to \infty$.

<!-- ======================= -->
<!-- PROBLEM 4.6             -->
<!-- ======================= -->
## Problem 4.6


The Lagrangian for two gravitating point masses is

$$
\mathcal{L} = \frac{1}{2}m_1 \dot{\mathbf{x}}_1^2 + \frac{1}{2}m_2 \dot{\mathbf{x}}_2^2 + \frac{Gm_1 m_2}{|\mathbf{x}_2 - \mathbf{x}_1|}.
$$

Define the center-of-mass and relative coordinates

$$
\mathbf{x}_{\mathrm{cm}} = \frac{m_1 \mathbf{x}_1 + m_2 \mathbf{x}_2}{M}, \qquad \mathbf{r} = \mathbf{x}_2 - \mathbf{x}_1,
$$

where $M = m_1 + m_2$. Inverting,

$$
\mathbf{x}_1 = \mathbf{x}_{\mathrm{cm}} - \frac{m_2}{M}\mathbf{r}, \qquad \mathbf{x}_2 = \mathbf{x}_{\mathrm{cm}} + \frac{m_1}{M}\mathbf{r}.
$$

The kinetic energy becomes

$$
\begin{align}
T &= \frac{1}{2}m_1\left(\dot{\mathbf{x}}_{\mathrm{cm}} - \frac{m_2}{M}\dot{\mathbf{r}}\right)^2 + \frac{1}{2}m_2\left(\dot{\mathbf{x}}_{\mathrm{cm}} + \frac{m_1}{M}\dot{\mathbf{r}}\right)^2 \\
&= \frac{1}{2}(m_1 + m_2)\dot{\mathbf{x}}_{\mathrm{cm}}^2 + \frac{1}{2}\frac{m_1 m_2}{M}\dot{\mathbf{r}}^2.
\end{align}
$$

Defining the reduced mass $\mu = m_1 m_2 / M$,

$$
\mathcal{L} = \frac{1}{2}M\dot{\mathbf{x}}_{\mathrm{cm}}^2 + \frac{1}{2}\mu\dot{\mathbf{r}}^2 + \frac{G\mu M}{|\mathbf{r}|}.
$$

The center-of-mass coordinate is cyclic ($\mathcal{L}$ does not depend on $\mathbf{x}_{\mathrm{cm}}$), so $M\dot{\mathbf{x}}_{\mathrm{cm}} = \text{const}$ and the first term is a constant of motion. The non-trivial dynamics is therefore governed by

$$
\mathcal{L}_{\mathrm{rel}} = \frac{1}{2}\mu\dot{\mathbf{r}}^2 + \frac{G\mu M}{|\mathbf{r}|},
$$

which is the Lagrangian of a particle of mass $\mu$ moving in the gravitational field of a fixed mass $M = m_1 + m_2$.

<!-- ======================= -->
<!-- PROBLEM 4.7             -->
<!-- ======================= -->
## Problem 4.7

We increase the mass with $M(t) = A(t) M_0$ where $A(t)$ is a smooth function that transitions from $1$ to $2$ around some time $t_{\rm mid}$, for example

$$
A(t) = \frac{3}{2} + \frac{1}{2}\tanh\!\left(\frac{t - t_{\rm mid}}{\tau}\right),
$$

so the potential at time $t$ is $\Phi(r,t) = A(t)\,\Phi_0(r)$ with $A$ going from $1$ to $2$. The parameter $\tau$ sets how fast the transition is compared to the orbital period $T_{\rm orb}$.

At each snapshot time $t_i$ we freeze the potential at its instantaneous amplitude $A(t_i)$ and compute the radial action

![Adiabatic Isochrone](assets/generated/p04_07_adiabatic_isochrone.png)

- **Left column (slow ramp, $\tau = 10\,T_{\rm orb}$):** The radial action $J_r$ stays nearly constant through the mass doubling — the adiabatic invariant is conserved. The orbit shrinks (smaller $R$) as the potential deepens, but does so in a way that preserves $J_r$.
- **Right column (fast ramp, $\tau = 0.1\,T_{\rm orb}$):** The mass doubles in a fraction of an orbital period. The orbit has no time to adjust adiabatically, and $J_r$ jumps to a new, larger value. The shaded band marks the ramp window.

<!-- ======================= -->
<!-- PROBLEM 4.8             -->
<!-- ======================= -->
## Problem 4.8

### Part a

Let's start from

$$
u'' + u = -\frac{F(1/u)}{L^2 u^2},
$$

where $u = 1/r$ and

$$
F(r) = -GM\left( \frac{1}{r^2} + \frac{2a}{r^3} \right).
$$

Therefore

$$
u'' + u = \frac{GM}{L^2}(1 + 2au),
$$

or equivalently

$$
u'' + \frac{1}{K^2}u = \frac{GM}{L^2} \qquad \text{where} \qquad \frac{1}{K^2} = 1 - \frac{2aGM}{L^2}. \tag{4.8.1}
$$

Solutions to this equation are of the form

$$
u = C \cos\left(\frac{\psi - \psi_0}{K}\right) + \frac{GM K^2}{L^2} = \frac{1}{r} \tag{4.8.2}
$$

To find the constant $C$, note that

$$
E = \frac{1}{2}\dot{r}^2 + \frac{L^2}{2r^2} + \Phi(r) = \frac{1}{2}\dot{r}^2 + \frac{L^2}{2r^2} - \frac{GM}{r} - \frac{aGM}{r^2}.
$$

We can write this in terms of $u$ by noting that

$$
\dot{r} = \frac{dr}{d\psi}\dot{\psi} = -\frac{1}{u^2}\frac{du}{d\psi}\frac{L}{r^2} = -L\frac{du}{d\psi},
$$

so that

$$
\begin{align}
L^2\left(\frac{du}{d\psi}\right)^2 &= 2E + 2GMu - (L^2 - 2aGM)u^2 \\
&\stackrel{(4.8.1)}{=} 2E + 2GMu - \frac{L^2}{K^2}u^2
\end{align}
$$

Using Eqn. (4.8.2) to substitute for $u$ and $du/d\psi$

$$
\begin{align}
\frac{L^2C^2}{K^2}\sin^2\left(\frac{\psi - \psi_0}{K}\right) &= 2E + 2GM\left(C \cos\left(\frac{\psi - \psi_0}{K}\right) + \frac{GM K^2}{L^2}\right) \\
&\quad - \frac{L^2}{K^2}\left(C \cos\left(\frac{\psi - \psi_0}{K}\right) + \frac{GM K^2}{L^2}\right)^2 \\
&= 2E + \frac{G^2M^2K^2}{L^2} - \frac{L^2 C^2}{K^2}\cos^2\left(\frac{\psi - \psi_0}{K}\right).
\end{align}
$$

Rearranging

$$
\frac{L^2C^2}{K^2} = 2E + \frac{G^2M^2K^2}{L^2} \quad \Rightarrow \quad E = \frac{L^2C^2}{2K^2} - \frac{1}{2}\left(\frac{GMK}{L}\right)^2
$$

### Part b

An integral of motion is any quantity $I(\mathbf{x}, \mathbf{v})$ that (i) depends only on the current phase-space coordinates and (ii) is conserved along the orbit. $\psi_0$ satisfies both conditions:

From Eqn. (4.8.2) and its derivative with respect to $\psi$,

$$
u - \frac{GMK^2}{L^2} = C\cos\left(\frac{\psi - \psi_0}{K}\right), \qquad \frac{du}{d\psi} = -\frac{C}{K}\sin\left(\frac{\psi - \psi_0}{K}\right),
$$

we can divide and solve for $\psi_0$:

$$
\psi_0 = \psi - K\mathrm{atan2}\left(-K\frac{du}{d\psi}, u - \frac{GMK^2}{L^2}\right).
$$

Since $u = 1/r$ and $du/d\psi = -\dot{r}/L$ (from Part a), this becomes

$$
\psi_0 = \psi - K\mathrm{atan2}\left(\frac{K\dot{r}}{L}, \dfrac{1}{r} - \dfrac{GMK^2}{L^2}\right),
$$

which expresses $\psi_0$ entirely in terms of the current phase-space coordinates $(r, \psi, \dot{r}, L)$, with no explicit time dependence. Since $\psi_0$ is also constant along the orbit by construction (it is a constant of integration of the orbit equation), it satisfies both conditions and is therefore an integral of motion.

This is analogous to the Laplace–Runge–Lenz vector in the Kepler problem: when $a = 0$ (so $K = 1$), $\psi_0$ gives the azimuthal angle of pericenter, encoding the same orientational information as the Runge–Lenz vector direction.

### Part c

For simplicity call $u_0 = GMK^2/L^2$, Eqn. (4.8.2) can be written as

$$
\cos\left(\frac{\psi - \psi_0}{K}\right) = \frac{1}{C}\left(u - u_0\right) \equiv A, \quad |A| \leq 1.
$$

Solutions are

$$
\psi = \psi_0 + K(\pm \arccos A + 2\pi n), \qquad n \in \mathbb{Z}. \tag{4.8.3}
$$

The key point here is that $\psi$ is an angle coordinate, so physically $\psi$ is only defined module $2\pi$, we then only care about solutions of Eqn. (4.8.3) mod $2\pi$. More explictly, call $\alpha = \arccos A \in [0, \pi]$, then question befomrs, how many distinct values are in the set

$$
S = \{ (\psi_0 + K(\pm \alpha + 2\pi n)) \pmod{2\pi} : n \in \mathbb{Z} \}.
$$

There are two branches here, $\psi^{\pm}_n = \psi_0 \pm K\alpha + 2\pi k n$, for each branch, increasing $n$ just adds the increment $2\pi K$, so the question above can be rephrased as: "On the circle $\pmod{2\pi}$, how many distinct points do we get by repeatedly adding $2\pi K$ to the initial point $\psi_0 \pm \alpha$?".

If $K \in \mathbb{Z}$, then $2\pi K n \pmod{2\pi} = 0$ for all $n$, and the entire family collapses to $\psi = \psi_0 \pm K \alpha \pmod{2\pi}$, that is, at most two distinct solutions physically this corresponds to the particle moving in, and moving out.

If, $K  = p/q\in \mathbb{Q}$, the are at most $2q$ distinct points.

Finally, if $K$ is irrational, the set is infinite. In summary we have

$$
|S| = \left\{
\begin{array}{ll}
\leq 2 & K \in \mathbb{Z} \\
\leq 2q & K = p/q \in \mathbb{Q} \\
\infty & K \not\in \mathbb{Q}
\end{array}
\right.
$$

So the condition "$(E, L, r, \psi_0)$ restricts $\psi$ to a finite set" forces $K \in \mathbb{Z}$.

### Part d

Inevetably we will only have two values of $\psi$ that satisfy the condition, but for $K = 1$ the radial period equal $2\pi$ and the orbit closes after on radial oscillation

$$
\frac{1}{K^2} = 1 = 1 - \frac{2aGM}{L^2} \quad \Rightarrow \quad a = 0,
$$

So we recover the Kepler potential.

<!-- ======================= -->
<!-- PROBLEM 4.10            -->
<!-- ======================= -->
## Problem 4.10 🌶️

Consider an orbit in an isochrone potential,

$$
\Phi(r) = -\frac{GM}{b + \sqrt{b^2 + r^2}},
$$

and suppose we apply an instantaneous velocity kick $\Delta \mathbf{v}$ at pericenter. We want to understand how the new pericenter compares to the original, and whether radial and tangential kicks differ in their effect.

### Part a: Orbit integration and pericenter measurement

```mermaid
flowchart TD
    A["Choose initial conditions (R, vR, vT)"] --> B["Integrate orbit in isochrone potential"]
    B --> C["Detect 1st pericenter via local minima in R(t)"]
    C --> D["Extract phase-space coords at pericenter"]
    D --> E["Pick next |dv| from sweep array"]
    E --> F{"Apply kick"}
    F -->|"Radial: vR -> vR + dv"| G["Integrate new orbit from kicked state"]
    F -->|"Tangential: vT -> vT + dv"| G
    G --> H["Measure new pericenter: min R(t)"]
    H --> I["Record ln(r_p_kicked / r_p)"]
    I --> J{"More |dv| values?"}
    J -->|Yes| E
    J -->|No| K["Plot results"]
```


The orbit comparison below shows the original orbit (solid) versus a radially kicked orbit (dashed), both starting from the same pericenter:

![Orbit Comparison](assets/generated/p04_10_orbit_comparison.png)

The plot below shows $\ln(r_p^{(k)} / r_p)$ versus kick magnitude, where $r_p$ is the original pericenter and $r_p^{(k)}$ is the pericenter after the kick:

![Parameter Study](assets/generated/p04_10_parameter_study.png)

- **Radial kicks** (left): the log-ratio is negative — the pericenter **decreases**. The orbit becomes more eccentric.
- **Tangential kicks** (right): the log-ratio is zero (within numerical noise).

### Part b: Pure tangential kick

From the equation

$$
\frac{1}{2}v_r^2 + \Phi_{\mathrm{eff}}(r) = E \quad \text{where} \quad \Phi_{\mathrm{eff}}(r;L_z) = \Phi(r) + \frac{L_z^2}{2r^2},
$$

we have that at the pericenter

$$
E = \Phi_{\mathrm{eff}}(r_p) = \Phi(r_p) + \frac{L_z^2}{2r_p^2}. \tag{4.10.1}
$$

Take $\Delta v_R = 0$ and $\Delta v_T > 0$, then immediately after the kick

$$
\begin{align}
L^{(k)} &= r_p (v_T + \Delta v_T) = L_z + r_p \Delta v_T, \\
E^{(k)} &= \frac{1}{2}(v_R^2 + (v_T + \Delta v_T)^2) + \Phi(r_p) = E + v_T \Delta v_T + \frac{1}{2}(\Delta v_T)^2. \tag{4.10.2}
\end{align}
$$

and the effective potential after the kick is

$$
\begin{align}
\Phi_{\mathrm{eff}}(r_p;L^{(k)}) &= \Phi(r_p) + \frac{(L^{(k)})^2}{2r_p^2} \\
&= \Phi(r_p) + \frac{L_z^2}{2r_p^2} + \frac{L_z}{r_p}\Delta v_T + \frac{(\Delta v_T)^2}{2} \\
&\stackrel{(4.10.1)}{=} E + v_T\Delta v_T + \frac{(\Delta v_T)^2}{2} \\
&\stackrel{(4.10.2)}{=} E^{(k)}.
\end{align}
$$

which means $r_p$ is still a turning point of the new orbit, so the pericenter does not change.

### Part b: Pure radial kick
Now take $\Delta v_T = 0$ and $\Delta v_R > 0$, then immediately after the kick

$$
\begin{align}
L^{(k)} &= r_p v_T = L_z, \tag{4.10.3} \\
E^{(k)} &= \frac{1}{2}((0 + \Delta v_R)^2 + v_T^2) + \Phi(r_p) = \Phi(r_p) + \frac{1}{2}v_T^2 + \frac{1}{2}(\Delta v_R)^2.
\end{align}
$$

Which means that

$$
\begin{align}
E^{(k)} - \Phi_{\mathrm{eff}}(r_p;L_z^{(k)}) &\stackrel{(4.10.3)}{=} E^{(k)} - \Phi_{\mathrm{eff}}(r_p;L_z) \\
&= \left[ \Phi(r_p) + \frac{1}{2}v_T^2 + \frac{1}{2}(\Delta v_R)^2 \right] - \left[ \Phi(r_p) + \frac{L_z^2}{2r_p^2} \right] \\
&= \frac{1}{2}(\Delta v_R)^2 > 0,
\end{align}
$$

So $r_p$ is no longer a turning point of the new orbit, the new pericenter is the new root of the equation $E^{(k)} = \Phi_{\mathrm{eff}}(r;L_z)$. To see in which direction the pericenter changes, note that $\Phi_{\mathrm{eff}}(r;L_z)$ is a decreasing function of $r$ for $r < r_p$ and an increasing function of $r$ for $r > r_p$. Since $E^{(k)} > \Phi_{\mathrm{eff}}(r_p;L_z)$, the new pericenter must be smaller than the original one, $r_p^{(k)} < r_p$, exactly as we found in the numerical experiment.

### Part c: Apocenter and eccentricity

The same analysis can be extended to the apocenter radius and the orbital eccentricity $e = (r_a - r_p)/(r_a + r_p)$:

![Apocenter Study](assets/generated/p04_10_apocenter_study.png)

![Eccentricity Study](assets/generated/p04_10_eccentricity_study.png)

- **Radial kicks** increase both the apocenter and the eccentricity — the orbit becomes larger and more elongated.
- **Tangential kicks** increase the apocenter but leave the eccentricity nearly unchanged, consistent with the pericenter being preserved.

!!! warning "Not sure where to start"
    I am not sure how to start this problem, I will come back to this later.
