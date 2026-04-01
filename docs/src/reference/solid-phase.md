
### Reference (2)

# Solid-Phase

To describe tidal and rotational deformations of a spherically symmetric body, `Obliqua` considers the spheroidal displacement–stress–gravity system. For each harmonic degree (``\ell``) and order (``m``), the spheroidal perturbed state is represented by the 6-vector

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

| Component       | Physical Meaning                                                                               | Units (SI) | Normalization Scale | Notes                           |
| --------------- | ---------------------------------------------------------------------------------------------- | ---------- | ------------------- | ------------------------------- |
| ``U_{\ell m}``    | Radial displacement                                                                            | m          | ``R_0``               | Typical planetary radius        |
| ``V_{\ell m}``    | Tangential displacement                                                                        | m          | ``R_0``               | Same as radial displacement     |
| ``R_{\ell m}``    | Radial stress                                                                                  | Pa         | ``\mu_0``             | Characteristic shear modulus    |
| ``S_{\ell m}``    | Tangential stress                                                                              | Pa         | ``\mu_0``             | Same as radial stress           |
| ``\Phi_{\ell m}`` | Gravitational potential perturbation                                                           | m``^2/``s``^2``      | ``g_0 R_0``           | ``g_0`` is characteristic gravity |
| ``Q_{\ell m}``    | Potential stress | m/s``^2``       | ``g_0``               | Same units as gravity           |

The spheroidal vector satisfies the first-order ODE system

```math
\frac{d\mathbf{y}_{\ell m}}{dr}
= \mathbf{A}_{\ell}(r)\,\mathbf{y}_{\ell m}(r).
```

Here, the coefficient matrix ``\mathbf{A}_{\ell}(r)`` represents the responds of the mantle to deformations, and is given by 

```math
\mathbf{A}_\ell(r) =
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

When solving the spheroidal displacement–stress–gravity system, the 6-vector can span vastly different physical units. This disparity can make the coefficient matrix ``\mathbf{A}_\ell(r)`` highly ill-conditioned, leading to numerical instability when computing linearly independent solutions or performing matrix inversions.

To mitigate this, we introduce a unit-normalization scaling matrix ``\mathbf{S}``, defined as
```math
  \mathbf{S} = \mathrm{diag}\Big(R_0, R_0, \mu_0, \mu_0, g_0 R_0, g_0 \Big),
```

where ``R_0, \mu_0, g_0`` are characteristic scales for length, stress, and gravity, respectively. The scaled variables are then

```math
  \tilde{\mathbf{y}}_{\ell m} = \mathbf{S}^{-1}\mathbf{y}_{\ell m}, \quad
  \tilde{\mathbf{A}}_\ell = \mathbf{S}^{-1} \mathbf{A}_\ell \mathbf{S}.
```

---

### Core–Mantle Boundary

In order to solve the system a general solution is constructed through propagation from the core-mantle boundary outwards. At the CMB radius ``r_C``, the spheroidal solution satisfies

```math
\mathbf{y}_\ell(r_C^+)
  = \mathbf{I}_C \mathbf{C},
```

where ``\mathbf{C} = (C_1, C_2, C_3)^T`` is a vector of integration constants determined by surface boundary conditions. The CMB Interface Matrix is given as 

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

The assumed normalization implies

```math
\mathbf{\tilde{I}}_C = \mathbf{S}^{-1} \mathbf{I}_C.
```

Once the constants ``\mathbf{C}`` are determined, the full perturbed state of the solid mantle is known.

---

### (Visco)elastic Solution

We propagate the solution using the so-called propagator matrix (``\pmb{\Pi}_\ell``). The propagator matrix solves the homogeneous differential system

```math
\frac{d\pmb{\Pi}_\ell(r, r')}{dr} = \pmb{A}_\ell(r)\,\pmb{\Pi}_\ell(r, r'),
```
at radius ``r`` w.r.t. the solution at the previous layer ``r'``, this is also know as the Cauchy data at radius (``r'``). Equivalently, we may write the normalized version of the propagator matrix as

```math
\frac{d\tilde{\pmb{\Pi}}_\ell(r, r')}{dr} = \tilde{\pmb{A}}_\ell(r)\,\tilde{\pmb{\Pi}}_\ell(r, r'),
```

Given that ``\tilde{\pmb{\Pi}}_\ell(r, r')`` is constructed from the normalized matrix ``\tilde{\pmb{A}}_\ell``, the normailization in ``\mathbf{A}_\ell`` enters ``\tilde{\pmb{\Pi}}_\ell`` as

```math
\tilde{\pmb{\Pi}}_\ell(r, r') = \left[\left[\left[  \mathbf{1} \times \mathbf{S}^{-1} \tilde{\pmb{\Pi}}_\ell(r_1, r_1') \mathbf{S}\right] \times \mathbf{S}^{-1} \tilde{\pmb{\Pi}}_\ell(r_2, r_2') \mathbf{S} \right] \times \dots \right] = \mathbf{S}^{-1} \pmb{\Pi}_\ell(r, r') \mathbf{S},
``` the normalization does not alter the physical solution. If ``r = r'`` we have

```math
\pmb{\Pi}_\ell(r', r') = \pmb{1}.
```

Each column of the propagator matrix is one of the six linearly independent solutions of

```math
\frac{d\pmb{y}_{\ell m}}{dr} = \pmb{A}_\ell(r)\,\pmb{y}_{\ell m}.
```

The six linearly independent solutions are multiplied by the propagator, forming a basis matrix. To prevent ill-conditioning due to widely differing units and growing/decaying solutions, we perform a QR decomposition at every sublayer:

```math
\pmb{\Pi}_\ell(r, r') = \mathbf{Q}\,\mathbf{R},
``` where ``\mathbf{Q}`` is an orthogonal matrix and ``\mathbf{R}`` is an upper triangular matrix. The orthogonal matrix ``\mathbf{Q}`` forms an orthonormalized basis used to propagate to the next sublayer, while the upper triangular matrix ``\mathbf{R}`` stores the mixing coefficients for back-propagation, ensuring accurate reconstruction of the physical solution. The process is repeated across all layers.

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

or in the normalized form

```math
\pmb{y}_{\ell m}(r) =
\mathbf{S} \, \tilde{\pmb{\Pi}}_\ell(r, r_C^+)\,\tilde{\pmb{I}}_C\,{\pmb{C}}.
```

This equation can be solved iteratively, up till the surface to yield the general responds of the interior to any form of tidal- or load induced deformation.

---

### Surface Boundary Conditions

Finally, to find the specific solution for the case where a planet is orbiting around a star-like object, we can specify the surface (``r=a^-``) boundary condition. In particular we need to constrain the 3rd, 4th, and 6th ``y``-values. To do this employ the projector on the 3rd, 4th, and 6th components given by 

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

Incase the solid part of the mantle extends up to the surface of the planet, then the surface is stress-free and the toroidal part has ``\sigma_{\ell m}^L = \pmb{0}``, and the solution no longer depends on ``\pmb{b}^L``. However, as in our formalism the solid is often loaded by a liquid magma ocean, we cannot ignore this term.

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

To solve this system we thus need only provide ``\pmb{P}_1\,\pmb{y}(a^-)``, we will provide some examples in Chapter 7 at the end of this component.

---

