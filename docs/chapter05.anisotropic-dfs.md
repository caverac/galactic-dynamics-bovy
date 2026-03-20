There are some expressions in section 5.6.4 that took me a while to figure out, I will write them here for my own reference and for anyone else who might find them useful.

<!-- ======================= -->
<!-- Eqn 5.79                -->
<!-- ======================= -->
## Deriving Eq. 5.79

To derive Eq. (5.79) note that the velocity can be written in spherical coordinates as $\mathbf{v} = (v_r, v_\theta, v_\phi)$, and we can further change $(v_\theta, v_\phi)$ to polar coordinates in the tangential velocity plane. Define the tangential speed as $v_t = \sqrt{v_\theta^2 + v_\phi^2}$, and the angle $\psi$ in the tangential plane so that $v_\theta = v_t \cos\psi$ and $v_\phi = v_t \sin\psi$. Then the volume element in velocity space is $d^3\mathbf{v} = (dv_\theta dv_\phi) dv_r  = (v_t\,dv_t\,d\psi)dv_r$.

$$
\begin{align}
\rho(r) &= \int d^3\mathbf{v} \; f(r, \mathbf{v}) = \int_0^{2\pi} d\psi \int dv_r dv_t\; v_t \, f(r, v_r, v_t) \\
&= 2\pi \int dv_r dv_t\; v_t \, f(r, v_r, v_t).
\end{align}
$$

Since $f = 0$ for $\mathcal{E} \geq 0$, we need to consider only $v_r^2 + v_t^2 \leq 2\Psi$. For a fixed $v_t$, $v_r$ ranges over $\pm\sqrt{2\Psi - v_t^2}$, and $v_t$ ranges from 0 to $\sqrt{2\Psi}$. Using the symmetry $f(\cdots, v_r) = f(\cdots, -v_r)$ to fold the $v_r$ integral:

$$
\rho(r) = 4\pi \int_0^{\sqrt{2\Psi}} dv_t \int_0^{\sqrt{2\Psi - v_t^2}} dv_r \; v_t \, f(\mathcal{E},\, r v_t)
$$

Now change the integration variable from $v_r$ to $\mathcal{E} = \Psi - \frac{1}{2}(v_r^2 + v_t^2)$
at fixed $v_t$. Then $d\mathcal{E} = -v_r\,dv_r$ and the limits are

- $v_r = 0 \;\Rightarrow\; \mathcal{E} = \Psi - v_t^2/2$
- $v_r = \sqrt{2\Psi - v_t^2} \;\Rightarrow\; \mathcal{E} = 0$

And the integral becomes

$$
\rho(r) = 4\pi \int_0^{\sqrt{2\Psi}} dv_t \int_0^{\Psi - v_t^2/2} d\mathcal{E} \; \frac{v_t \, f(\mathcal{E},\, r v_t)}{\sqrt{2(\Psi - \mathcal{E}) - v_t^2}}
$$

Finally, swap the order of integration. The region is
$0 \leq v_t \leq \sqrt{2\Psi}$ and $0 \leq \mathcal{E} \leq \Psi - v_t^2/2$,
which is equivalent to $0 \leq \mathcal{E} \leq \Psi$ and
$0 \leq v_t \leq \sqrt{2(\Psi - \mathcal{E})}$:

$$
\boxed{
\rho(r) = 4\pi \int_0^{\Psi} d\mathcal{E} \int_0^{\sqrt{2(\Psi - \mathcal{E})}} dv_t \; \frac{v_t \, f(\mathcal{E},\, r v_t)} {\sqrt{2(\Psi - \mathcal{E}) - v_t^2}}
}
\tag{5.79}
$$

<!-- ======================= -->
<!-- Eqn 5.81                -->
<!-- ======================= -->
## Deriving Eq. 5.81 from Eq. 5.79

To obtain $\rho\,\sigma_r^2$ we need the second radial-velocity moment:

$$
\rho(r)\,\sigma_r^2(r) = \int d^3\mathbf{v}\; v_r^2 f(r, \mathbf{v})
$$

Repeating the same steps as above (integrate $\psi$, fold $v_r$, change variable to $\mathcal{E}$), the only difference is the extra $v_r^2$ in the integrand. Since $v_r^2\,dv_r = v_r\,d\mathcal{E}$, the $1/v_r$ Jacobian and $v_r^2$ moment combine into a single $v_r = \sqrt{2(\Psi-\mathcal{E})-v_t^2}$ in the **numerator** (instead of the denominator as in 5.79). After swapping integration order:

$$
\boxed{
\rho(r)\,\sigma_r^2(r) = 4\pi \int_0^{\Psi} d\mathcal{E}
    \int_0^{\sqrt{2(\Psi - \mathcal{E})}} dv_t \;
    v_t \, f(\mathcal{E},\, r v_t)\,
    \sqrt{2(\Psi - \mathcal{E}) - v_t^2}
}
\tag{5.81}
$$

<!-- ======================= -->
<!-- Eqn 5.82                -->
<!-- ======================= -->
## Deriving Eq. 5.82 from Eq. 5.81

Equation 5.80 defines the augmented density as a function of $(\Psi, r)$ - identical in form to the density integral (5.79), but treating $\Psi$ as a free variable rather than evaluating it at $\Psi[r]$:

$$
\tilde{\rho}(\Psi, r) = 4\pi \int_0^{\Psi} d\mathcal{E} \int_0^{\sqrt{2(\Psi - \mathcal{E})}} dv_t \; \frac{v_t \, f(\mathcal{E},\, r v_t)} {\sqrt{2(\Psi - \mathcal{E}) - v_t^2}} \tag{5.80}
$$

with the true density recovered as $\rho(r) = \tilde{\rho}(\Psi[r], r)$.

Now consider Eq. 5.81

$$
\begin{align}
\tilde{\rho}(\Psi, r)\,\sigma_r^2(r) &= 4\pi \int_0^{\Psi} d\mathcal{E} \int_0^{\sqrt{2(\Psi-\mathcal{E})}} dv_t \; v_t\, f(\mathcal{E}, rv_t)\, \sqrt{2(\Psi-\mathcal{E})-v_t^2} \\
&= 4\pi \int_0^{\Psi} d\mathcal{E} \int_0^{\sqrt{2(\Psi-\mathcal{E})}} dv_t \; v_t\, f(\mathcal{E}, rv_t)\, \int_{\mathcal{E} + v_t^2/2}^{\Psi} \frac{d\Psi'}{\sqrt{2(\Psi'-\mathcal{E})-v_t^2}} \\
&= 4\pi \int_0^{\Psi} d\Psi' \int_0^{\Psi'} d\mathcal{E}\int_0^{\sqrt{2(\Psi'-\mathcal{E})}} dv_t \; \frac{v_t\, f(\mathcal{E}, rv_t)}{\sqrt{2(\Psi'-\mathcal{E})-v_t^2}} \\
&\stackrel{5.80}{=} \int_0^{\Psi} d\Psi' \; \tilde{\rho}(\Psi', r)
\end{align}
$$

Where in the secont-to-last step we have re-arranged the domain

$$
\begin{align}
0 \leq \mathcal{E} \leq \Psi&, \quad 0 \leq v_t \leq \sqrt{2(\Psi-\mathcal{E})}, \quad \mathcal{E} + \frac{v_t^2}{2} \leq \Psi' \leq \Psi \\
\Rightarrow\quad& 0 \leq \Psi' \leq \Psi, \quad 0 \leq \mathcal{E} \leq \Psi', \quad 0 \leq v_t \leq \sqrt{2(\Psi'-\mathcal{E})}.
\end{align}
$$

Thus we have

$$
\tilde{\rho}(\Psi, r) = \frac{1}{\sigma_r^2(r)} \frac{d}{d\Psi}\tilde{\rho}(\Psi, r) \tag{5.82}
$$
