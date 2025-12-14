
### Reference (2)

# Chapter 2 - Solid-Phase

To describe tidal and rotational deformations of a spherically symmetric body, `fwlLove` considers the spheroidal displacement–stress–gravity system. For each harmonic degree (``\ell``) and order (``m``), the spheroidal perturbed state is represented by the 6-vector

```math
    \mathbf{y}_{\ell m} =
    \begin{pmatrix}
    U_{\ell m} \
    V_{\ell m} \
    R_{\ell m} \
    S_{\ell m} \
    \Phi_{\ell m} \
    Q_{\ell m}
    \end{pmatrix},
```

where

* ``U_{\ell m}``: radial displacement
* ``V_{\ell m}``: tangential displacement
* ``R_{\ell m}``: radial stress
* ``S_{\ell m}``: tangential stress
* ``\Phi_{\ell m}``: gravitational potential perturbation
* ``Q_{\ell m}``: “potential stress,” defined as

```math
Q_{\ell m}
= \frac{\partial \Phi_{\ell m}}{\partial r}
  + \frac{\ell+1}{r}\Phi_{\ell m}
  + 4\pi G \rho_0 U_{\ell m}.
```

### Governing System

The spheroidal vector satisfies the first-order ODE system

```math
\frac{d\mathbf{y}_{\ell m}}{dr}
= \mathbf{A}_{\ell}(r)\,\mathbf{y}_{\ell m}(r).
```

The coefficient matrix ``\mathbf{A}_{\ell}(r)`` is

```math
A(r) =
\begin{pmatrix}
-\frac{2\lambda}{r\beta} &
\frac{\ell(\ell+1)\lambda}{r\beta} &
\frac{1}{\beta} & 0 & 0 & 0 \\[1.2em]

-\frac{1}{r} &
\frac{1}{r} &
0 &
\frac{1}{\mu} &
0 & 0 \\[1.2em]

\frac{4}{r}\!\left( \frac{3\kappa\mu}{r\beta} - \rho_{0}g \right) - \rho_{0}\omega^{2} &
\frac{\ell(\ell+1)}{r}\!\left(\rho_{0}g - \frac{6\kappa\mu}{r\beta}\right) &
-\frac{4\mu}{r\beta} &
\frac{\ell(\ell+1)}{r} &
-\frac{\rho_{0}(\ell+1)}{r} &
\rho_{0} \\[1.2em]

\frac{1}{r}\!\left(\rho_{0}g - \frac{6\mu\kappa}{r\beta}\right) &
\frac{2\mu}{r^{2}}\!\left[\ell(\ell+1)\!\left(1+\frac{\lambda}{\beta}\right)-1\right] - \rho_{0}\omega^{2} &
-\frac{\lambda}{r\beta} &
-\frac{3}{r} &
\frac{\rho_{0}}{r} &
0 \\[1.2em]

-4\pi G \rho_{0} &
0 & 0 & 0 &
-\frac{\ell+1}{r} &
1 \\[1.2em]

-\frac{4\pi G \rho_{0} (\ell+1)}{r} &
\frac{4\pi G \rho_{0} \ell(\ell+1)}{r} &
0 & 0 & 0 &
\frac{\ell-1}{r}
\end{pmatrix}.
```

The material parameters satisfy

```math
\beta = \lambda + 2\mu
```

and

```math
\lambda = \kappa - \frac{2}{3}\mu.
```

The background fields ``\rho_0``, ``g`` and ``\kappa`` are obtained from SPIDER/ARAGOG.

---

### Core–Mantle Boundary (CMB)

At the CMB radius ``r_C``, the spheroidal solution satisfies

```math
\mathbf{y}_\ell(r_C^+)
  = \mathbf{I}_C \mathbf{C},
```

where ``\mathbf{C} = (C_1, C_2, C_3)^T`` is a vector of integration constants determined by surface boundary conditions.

### CMB Interface Matrix

```math
\mathbf{I}_C =
\begin{pmatrix}
-\psi_\ell(r_C)/g(r_C) & 0 & 1 \\[1.2em]
0 & 1 & 0 \\[1.2em]
0 & 0 & g(r_C)/\rho_0(r_C^-) \\[1.2em]
0 & 0 & 0 \\[1.2em]
\psi_\ell(r_C) & 0 & 0 \\[1.2em]
q_\ell(r_C) & 0 & 4\pi G \rho_0(r_C^-)
\end{pmatrix}.
```

Once the constants ``\mathbf{C}`` are determined, the full perturbed state of the solid mantle is known.

---

### (Visco)elastic Solution

The propagator matrix (``\pmb{\Pi}_\ell``) (a ``6\times 6`` matrix for the spheroidal case) solves the homogeneous differential system

```math
\frac{d\pmb{\Pi}_\ell(r, r')}{dr} = \pmb{A}_\ell(r)\,\pmb{\Pi}_\ell(r, r'),
```

with Cauchy data at radius (``r'``):

```math
\pmb{\Pi}_\ell(r', r') = \pmb{1}.
```

Each column of the propagator matrix is one of the six linearly independent solutions of

```math
\frac{d\pmb{y}_{\ell m}}{dr} = \pmb{A}_\ell(r)\,\pmb{y}_{\ell m}.
```

We impose continuity:

```math
\pmb{\Pi}_\ell(r_j^+, r') = \pmb{\Pi}_\ell(r_j^-, r'),
```

and apply CMB boundary conditions:

```math
\pmb{y}_{\ell m}(r_C^+) = \pmb{y}_0 = \pmb{I}_C\,\pmb{C}.
```

Therefore,

```math
\pmb{y}_{\ell m}(r) =
\pmb{\Pi}_\ell(r, r_C^+)\,\pmb{I}_C\,\pmb{C}.
```

---

### Surface Boundary Conditions

The projector on the 3rd, 4th, and 6th components is

```math
\pmb{P}_1\,\pmb{y}(a^-) =
\pmb{P}_1 \left[
\pmb{\Pi}_\ell(a^-, r_C^+)\,\pmb{I}_C\,\pmb{C}
\right]
= \pmb{b},
```

where

```math
\pmb{b} = \sigma_{\ell m}^L\,\pmb{b}^L +
\left(\Phi_{\ell m}^T(a) + \Phi_{\ell m}^C(a)\right)\pmb{b}^T,
```

with

```math
\pmb{b}^L =
\begin{pmatrix}
-\dfrac{(2\ell+1)g(a)}{4\pi a^2} \\[1em]
0 \\[1em]
-\dfrac{(2\ell+1)G}{a^2}
\end{pmatrix}
\qquad\text{(Load)},
```

```math
\pmb{b}^T =
\begin{pmatrix}
0 \\[1em]
0 \\[1em]
\dfrac{2\ell+1}{a}
\end{pmatrix}
\qquad\text{(Tidal)}.
```

The toroidal part has (``\pmb{b} = \pmb{0}``) for a stress-free surface.

The integration constants follow from

```math
\pmb{C} =
\left(\pmb{P}_1 \pmb{\Pi}_\ell(a, r_C)\,\pmb{I}_C\right)^{-1}\pmb{b}.
```

Thus,

```math
\pmb{y}_{\ell m}(r)
=
\pmb{\Pi}_\ell(r, r_C)\,
\pmb{I}_C\,
\left(\pmb{P}_1 \pmb{\Pi}_\ell(a, r_C)\,\pmb{I}_C\right)^{-1}
\pmb{b}.
```

---

### Semidiurnal Lunar Tide

For (``\ell = 2``), the semidiurnal case:

* The boundary conditions at the surface are
  (``y_3 = y_4 = 0``)
  (see *Takeuchi et al. 1972*).

The final condition is

```math
y_6 =
- \frac{(2\ell + 1) g_s}{4}
\frac{M_M}{M_E}
\left(\frac{R_E}{a}\right)^3.
```

Thus,

```math
\pmb{P}_1\,\pmb{y}(a^-) =
\begin{pmatrix}
0 \\[1em]
0 \\[1em]
-\dfrac{(2\ell+1) g_s}{4}
\dfrac{M_M}{M_E}
\left(\dfrac{R_E}{a}\right)^3
\end{pmatrix}.
```
---


