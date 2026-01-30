# Chapter 3: Elements of Classical Mechanics

<!-- ======================= -->
<!-- PROBLEM 3.1             -->
<!-- ======================= -->
## Problem 3.1

If the mass of the Milky Way $M$ is concentrated within the Sun's orbit $R_0$, then applying Newton's first theorem we have

$$
\Phi(r) = -\frac{GM}{r}
$$

and the circular velocity is

$$
v_c(r) = \sqrt{r\frac{d\Phi}{dr}} = \sqrt{\frac{GM}{r}}
$$

The scape velocity is

$$
v_{\rm esc}(R_0) = \sqrt{-2\Phi(R_0)} = \sqrt{\frac{2GM}{R_0}} = \sqrt{2}v_c(R_0)
$$

<!-- ======================= -->
<!-- PROBLEM 3.2             -->
<!-- ======================= -->
## Problem 3.2

At a very high-level, we have an optimzation problem with

- Contraints: $M(< 8\, \mathrm{kpc})$, $v_{\rm esc}(8\, \mathrm{kpc})$, and the Dutton & Maccio (2014) relation $M_{200} = M_{200}(c)$.
- Variables: $\rho_0$ and $a$.

That is over-constrained! What it means for practical purpose is that we can only satisfy two of the three constraints at the same time.

### Using Dutton & Maccio (2014)
For an NFW potential the mass is

$$
M(<r) = M_{200}\frac{f(r/a)}{f(c)}. \tag{3.2.1}
$$

From the estimation $M(r < 8\, \mathrm{kpc}) \approx 9\times 10^{10} M_\odot$ we have

$$
9\times 10^{10} M_\odot = M_{200}(c)\frac{f(8\, \mathrm{kpc} \cdot c/R_{200})}{f(c)}, \tag{3.2.2}
$$

where I have explicitly indicated that $M_{200}$ depends on the concentration $c$ through the relation

$$
\log c = 0.905 - 0.101\log\left(\frac{M_{200}}{10^{12}h^{-1}M_\odot}\right).
$$

Solving Eqn. (3.2.2) for $c$ we get (note that $R_{200}$ also depends on $M_{200}$ and therefore on $c$)

```python
>>> from galactic_dynamics_bovy.chapter03.nfw_concentration import solve_concentration
>>> solve_concentration(m_inner=9e10, r_inner=8.0)
SolvedNFWModel(c=5.1647268050004875, m200=11869.883211996499, r200=1043.8048483715481, a=202.102625711187, v_esc_predicted=1799.0994645175906)
```

This is clearly a horrible prediction for the escape velocity at the Sun's position!

$$
M_{200} \approx 1.2 \times 10^{14} M_\odot, \quad c \approx 5.16
$$

### Ignoring Dutton & Maccio (2014)
If we ignore the Dutton & Maccio (2014) relation and instead use the escape we can proceed as in the example in the notes. We will have two equations, and two unknowns: $M_{200}$ and $c$. The constraints are $M(< 8\, \mathrm{kpc}) = 9\times 10^{10} M_\odot$ and $v_{\rm esc}(8\, \mathrm{kpc}) = 550\, \mathrm{km/s}$.


For an NFW halo truncated at $R_{200}$, the potential at radius $r$ has two contributions:

$$
\Phi(r) = \underbrace{-\frac{GM(<r)}{r}}_{\text{mass inside } r} \underbrace{- \frac{GM_{200}}{f(c)}\left[\frac{1}{a+r} - \frac{1}{a+R_{200}}\right]}_{\text{shell from } r \text{ to } R_{200}} \tag{3.2.3}
$$

where $a = R_{200}/c$ is the scale radius and $f(x) = \ln(1+x) - x/(1+x)$.

At $R_{200}$, the second term vanishes and we simply have $\Phi(R_{200}) = -GM_{200}/R_{200}$.

The escape velocity from $r$ to $R_{200}$ is then

$$
v_{\rm esc}^2(r) = 2\left[\Phi(R_{200}) - \Phi(r)\right] = \frac{2GM(<r)}{r} + \frac{2GM_{200}}{f(c)}\left[\frac{1}{a+r} - \frac{1}{a+R_{200}}\right] - \frac{2GM_{200}}{R_{200}} \tag{3.2.4}
$$

This gives us two equations for two unknowns $(c, M_{200})$:

$$
\begin{cases}
M(<r) = M_{200}\dfrac{f(r/a)}{f(c)} = 9\times 10^{10} M_\odot \\[1em]
v_{\rm esc}(r) = 550\, \mathrm{km/s}
\end{cases}
$$

Solving these two equations simultaneously we get

```python
>>> from galactic_dynamics_bovy.chapter03.nfw_concentration import solve_concentration_from_vesc
>>> solve_concentration_from_vesc(m_inner=9e10, r_inner=8.0, v_esc=550.0)
SolvedNFWModel(c=25.7494572899753, m200=123.74644450054568, r200=228.02444341160776, a=8.855504830402058, v_esc_predicted=550.0000000000005)
```

This is more reasonable

$$
M_{200} \approx 1.2 \times 10^{12} M_\odot, \quad c \approx 25.7
$$

however, the concentration does not agree with the Dutton & Maccio (2014) relation.

!!! question "Discussion"
    I was not able to find a good compromise between the two constraints and the Dutton & Maccio (2014) relation. Maybe I am doing something extraordinarily wrong? If you find a better solution please let me know!


<!-- ======================= -->
<!-- PROBLEM 3.3             -->
<!-- ======================= -->
## Problem 3.3

In spherical coordinates, the velocity vector can be written as

$$
\mathbf{v} = \dot{r}\hat{r} + r\dot{\theta}\hat{\theta} + r\sin\theta\dot{\phi}\hat{\phi},
$$

so that the kinetic energy is

$$
T = \frac{1}{2}m\left(\dot{r}^2 + r^2\dot{\theta}^2 + r^2\sin^2\theta\dot{\phi}^2\right).
$$

In a conservative potential $\partial \mathcal{L} /\partial p_j = \partial \mathcal{T} /\partial p_j$, therefore the generalized momenta are

$$
\begin{align}
\frac{\partial \mathcal{L}}{\partial \dot{r}} &= m\dot{r} = p_r, \\
\frac{\partial \mathcal{L}}{\partial \dot{\theta}} &= mr^2\dot{\theta} = p_\theta, \\
\frac{\partial \mathcal{L}}{\partial \dot{\phi}} &= mr^2\sin^2\theta\dot{\phi} = p_\phi.
\end{align}
$$

<!-- ======================= -->
<!-- PROBLEM 3.4             -->
<!-- ======================= -->
## Problem 3.4

For $F = F_2(\mathbf{q}, \mathbf{p'}, t) - \mathbf{q}'\cdot\mathbf{p'}$ the equation

$$
\dot{\mathbf{q}}\cdot\mathbf{p} - H = \dot{\mathbf{q}'}\cdot\mathbf{p'} - K + \frac{dF}{dt}
$$

becomes

$$
\dot{\mathbf{q}}\cdot\mathbf{p} - H = \dot{\mathbf{q}'}\cdot\mathbf{p'} - K + \frac{\partial F_2}{\partial t} + \dot{\mathbf{q}}\cdot\frac{\partial F_2}{\partial \mathbf{q}} + \dot{\mathbf{p}'}\cdot\frac{\partial F_2}{\partial \mathbf{p'}} - \dot{\mathbf{q}'}\cdot\mathbf{p'} - \mathbf{q}'\cdot\dot{\mathbf{p}'}.
$$

Rearranging to group terms by their time derivatives:

$$
\dot{\mathbf{q}}\cdot\left(\mathbf{p} - \frac{\partial F_2}{\partial \mathbf{q}}\right) + \dot{\mathbf{p}'}\cdot\left(\mathbf{q}' - \frac{\partial F_2}{\partial \mathbf{p}'}\right) + \left(K - H - \frac{\partial F_2}{\partial t}\right) = 0.
$$

Since $\dot{\mathbf{q}}$, $\dot{\mathbf{p}'}$, and the time evolution are independent, each coefficient must vanish:

$$
\mathbf{p} = \frac{\partial F_2}{\partial \mathbf{q}}, \quad \mathbf{q}' = \frac{\partial F_2}{\partial \mathbf{p}'}, \quad K = H + \frac{\partial F_2}{\partial t}.
$$


<!-- ======================= -->
<!-- PROBLEM 3.5             -->
<!-- ======================= -->
## Problem 3.5

### Part a
For $F_2(\mathbf{q}, \mathbf{p'}, t) = \mathbf{q}\cdot\mathbf{p'}$ we have

$$
\begin{align}
\mathbf{p} &= \frac{\partial F_2}{\partial \mathbf{q}} = \mathbf{p'}, \\
\mathbf{q}' &= \frac{\partial F_2}{\partial \mathbf{p'}} = \mathbf{q}.
\end{align}
$$

$(\mathbf{q}, \mathbf{p}) \to (\mathbf{q}', \mathbf{p}') = (\mathbf{q}, \mathbf{p})$

### Part b
For $F_2(\mathbf{q}, \mathbf{p'}, t) = \mathbf{f}(\mathbf{q})\cdot\mathbf{p'}$ we have

$$
\begin{align}
\mathbf{p} &= \frac{\partial F_2}{\partial \mathbf{q}} = \frac{\partial \mathbf{f}}{\partial \mathbf{q}}\cdot\mathbf{p'}, \\
\mathbf{q}' &= \frac{\partial F_2}{\partial \mathbf{p'}} = \mathbf{f}(\mathbf{q}).
\end{align}
$$

$(\mathbf{q}, \mathbf{p}) \to (\mathbf{q}', \mathbf{p}') = (\mathbf{f}(\mathbf{q}), [\partial \mathbf{f}/\partial \mathbf{q}]^{-1}\cdot\mathbf{p})$, provided the inverse exists.

### Part c

For $F = F_3(\mathbf{q}', \mathbf{p}, t) + \mathbf{q}\cdot\mathbf{p}$ the equation

$$
\dot{\mathbf{q}}\cdot\mathbf{p} - H = \dot{\mathbf{q}'}\cdot\mathbf{p'} - K + \frac{dF}{dt}
$$

becomes

$$
\dot{\mathbf{q}}\cdot\mathbf{p} - H = \dot{\mathbf{q}'}\cdot\mathbf{p'} - K + \frac{\partial F_3}{\partial t} + \dot{\mathbf{q}'}\cdot\frac{\partial F_3}{\partial \mathbf{q}'} + \dot{\mathbf{p}}\cdot\frac{\partial F_3}{\partial \mathbf{p}} + \dot{\mathbf{q}}\cdot\mathbf{p} + \mathbf{q}\cdot\dot{\mathbf{p}}.
$$

Rearranging a bit

$$
\dot{\mathbf{q}'}\cdot\left(\mathbf{p'} + \frac{\partial F_3}{\partial \mathbf{q}'}\right) + \dot{\mathbf{p}}\cdot\left(\mathbf{q} + \frac{\partial F_3}{\partial \mathbf{p}}\right) + \left(K - H - \frac{\partial F_3}{\partial t}\right) = 0.
$$

Since $\dot{\mathbf{q}'}$, $\dot{\mathbf{p}}$, and the time evolution are independent, each coefficient must vanish:

$$
\mathbf{p'} = -\frac{\partial F_3}{\partial \mathbf{q}'}, \quad \mathbf{q} = -\frac{\partial F_3}{\partial \mathbf{p}}, \quad K = H + \frac{\partial F_3}{\partial t}.
$$

Now take $F_3(\mathbf{q}', \mathbf{p}, t) = \mathbf{f}(\mathbf{q}')\cdot\mathbf{p}$

$$
\begin{align}
\mathbf{p'} &= -\frac{\partial F_3}{\partial \mathbf{q}'} = -\frac{\partial \mathbf{f}}{\partial \mathbf{q}'}\cdot\mathbf{p}, \\
\mathbf{q} &= -\frac{\partial F_3}{\partial \mathbf{p}} = -\mathbf{f}(\mathbf{q}').
\end{align}
$$

<!-- ======================= -->
<!-- PROBLEM 3.6             -->
<!-- ======================= -->
## Problem 3.6

### Hamilton-Jacobi Equation

The equation we need to solve is

$$
\left(\frac{d W}{d x}\right)^2 + \omega^2 x^2 = 2E,
$$

which can be integrated once

$$
\begin{align}
W(x; E) &= \pm \int dx \sqrt{2E - \omega^2 x^2} \quad \cos \psi = -\frac{\omega x}{\sqrt{2E}} \\
&= \pm \int d\psi\, \left(\frac{\sqrt{2E}}{\omega}\sin \psi\right) (\sqrt{2E}\sin \psi) \\
&= \pm \frac{2E}{\omega} \int d\psi\,\sin^2 \psi \\
&= \pm \frac{E}{2\omega} \left(2\psi - \sin 2\psi\right) + C
\end{align}
$$

To calculate the action $J$ we need to note that over a full orbit $x$ goes from $-x_{\rm max}$ to $+x_{\rm max}$ and back, where $x_{\rm max} = \sqrt{2E}/\omega$

$$
\begin{align}
J &= \frac{1}{2\pi} \oint dx\, p \\
&= \frac{1}{2\pi} 2\int_{-x_{\rm max}}^{x_{\rm max}} dx\, \sqrt{2E - \omega^2 x^2} \\
&= \frac{E}{\omega \pi} \int_0^{2\pi} d\psi\, \sin^2 \psi \\
&= \frac{E}{\omega}
\end{align}
$$

To calculate the frequency note that on the orbit $H(J) = E(J) = \omega J$, therefore

$$
\Omega = \frac{dH(J)}{dJ} = \omega.
$$

This is the expected result for a harmonic oscillator! Also, note that

$$
W(J, \psi) = J\left(\psi - \frac{1}{2}\sin 2\psi\right) + C \qquad\text{with}\qquad \cos \psi = -x\sqrt{\frac{\omega}{2J}}.
$$

Therefore,

$$
\begin{align}
\theta = \left.\frac{\partial W}{\partial J}\right|_x &= \frac{\partial W}{\partial J} + \left.\frac{\partial W}{\partial \psi}\frac{\partial \psi}{\partial J} \right|_x \\
&= \left(\psi - \frac{1}{2}\sin 2\psi\right) + J (1 - \cos 2\psi)\left(\frac{1}{2J}\cot \psi\right) \\
&= \psi -\frac{1}{2}\sin 2\psi + \frac{1}{2}(2\sin^2\psi)\frac{\cos\psi}{\sin\psi} \\
&= \psi - \frac{1}{2}\sin 2\psi + \frac{1}{2}\sin 2\psi \\
&= \psi.
\end{align}
$$

That is $\psi$ is indeed the angle variable!

We can write now $x$ and $p$ in terms of the action-angle variables:

$$
\begin{align}
x &= -\sqrt{\frac{2J}{\omega}}\cos \theta, \\
p &= \frac{dW}{dx} = \sqrt{2E - \omega^2x^2} = \sqrt{2J\omega}\sin \theta.
\end{align}
$$

### Multidimensional case

For a 2D harmonic oscillator with separable potential

$$
H = \frac{1}{2}p_x^2 + \frac{1}{2}p_y^2 + \frac{1}{2}\omega_x^2 x^2 + \frac{1}{2}\omega_y^2 y^2,
$$

the Hamilton-Jacobi equation is

$$
\frac{1}{2}\left(\frac{\partial W}{\partial x}\right)^2 + \frac{1}{2}\left(\frac{\partial W}{\partial y}\right)^2 + \frac{1}{2}\omega_x^2 x^2 + \frac{1}{2}\omega_y^2 y^2 = E.
$$

Since the potential is separable, we try the ansatz $W(x, y) = W_x(x) + W_y(y)$:

$$
\underbrace{\frac{1}{2}\left(\frac{dW_x}{dx}\right)^2 + \frac{1}{2}\omega_x^2 x^2}_{\text{function of } x \text{ only}} + \underbrace{\frac{1}{2}\left(\frac{dW_y}{dy}\right)^2 + \frac{1}{2}\omega_y^2 y^2}_{\text{function of } y \text{ only}} = E.
$$

Rearranging:

$$
\frac{1}{2}\left(\frac{dW_x}{dx}\right)^2 + \frac{1}{2}\omega_x^2 x^2 = E - \frac{1}{2}\left(\frac{dW_y}{dy}\right)^2 - \frac{1}{2}\omega_y^2 y^2.
$$

The LHS depends only on $x$, the RHS depends only on $y$. For this equality to hold for all $x$ and $y$, both sides must equal a constant $E_x$:

$$
\begin{cases}
\dfrac{1}{2}\left(\dfrac{dW_x}{dx}\right)^2 + \dfrac{1}{2}\omega_x^2 x^2 = E_x \\[1em]
\dfrac{1}{2}\left(\dfrac{dW_y}{dy}\right)^2 + \dfrac{1}{2}\omega_y^2 y^2 = E_y
\end{cases}
\quad \text{with} \quad E_x + E_y = E.
$$

Each equation is identical to the 1D case solved above! Therefore:

$$
W(x, y) = W_x(x; E_x) + W_y(y; E_y),
$$

with actions $J_x = E_x/\omega_x$ and $J_y = E_y/\omega_y$, and frequencies $\Omega_x = \omega_x$ and $\Omega_y = \omega_y$.

The phase-space coordinates become:

$$
\begin{align}
x &= -\sqrt{\frac{2J_x}{\omega_x}}\cos \theta_x, & p_x &= \sqrt{2J_x\omega_x}\sin \theta_x, \\
y &= -\sqrt{\frac{2J_y}{\omega_y}}\cos \theta_y, & p_y &= \sqrt{2J_y\omega_y}\sin \theta_y.
\end{align}
$$

This separation works for any potential of the form $\Phi(x, y) = \Phi_x(x) + \Phi_y(y)$, not just harmonic oscillators.
