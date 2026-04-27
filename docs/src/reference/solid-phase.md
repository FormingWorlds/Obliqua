
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

Note: the density ``\rho_0`` used here corresponds to the core density. The assumed normalization implies

```math
\mathbf{\tilde{I}}_C = \mathbf{S}^{-1} \mathbf{I}_C.
```

Once the constants ``\mathbf{C}`` are determined, the full perturbed state of the solid mantle is known.

---

### (Visco)elastic Solution

We provide two methods to solve the coupled system of ODEs, the shooting method and the relaxation method.

#### Shooting method

propagate the solution using the so-called propagator matrix (``\pmb{\Pi}_\ell``). The propagator matrix solves the homogeneous differential system

```math
\frac{d\pmb{\Pi}_\ell(r, r')}{dr} = \pmb{A}_\ell(r)\,\pmb{\Pi}_\ell(r, r'),
```
at radius ``r`` w.r.t. the solution at the previous layer ``r'``, this is also know as the Cauchy data at radius (``r'``). If ``r = r'`` we have

```math
\pmb{\Pi}_\ell(r', r') = \pmb{1}.
```

Each column of the propagator matrix is one of the six linearly independent solutions of

```math
\frac{d\pmb{y}_{\ell m}}{dr} = \pmb{A}_\ell(r)\,\pmb{y}_{\ell m}.
```

The six linearly independent solutions are multiplied by the propagator, forming a basis matrix. We impose continuity:

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

This equation can be solved iteratively, up till the surface to yield the general responds of the interior to any form of tidal- or load induced deformation.

#### Relaxation method

Alternatively, to solve the above equation numerically, we approximate it using a second-order finite-difference scheme:

```math
\mathbf{y}_{n+1} - \mathbf{y}_n = \frac{\Delta r}{2} \left( \mathbf{A}_{n+1} \mathbf{y}_{n+1} + \mathbf{A}_n \mathbf{y}_n \right)
```

Here, ``\mathbf{y}_n`` and ``\mathbf{A}_n`` are evaluated at radius ``r = r_n`` for ``n = 1, \dots, N``, and ``\Delta r = r_{n+1} - r_n`` is the grid spacing.

The grid is non-uniform and follows an exponential distribution, with the largest spacing ``\Delta r_{\max}`` located just above the core–mantle boundary (CMB), and the smallest spacing ``\Delta r_{\min}`` near the surface.

The total number of nodes is given by:

```math
N = \mathrm{round} \left(
\frac{R_E - R_c}{\Delta r_{\min}} \cdot \frac{\alpha}{\exp(\alpha) - 1}
\right)
```

where ``R_c`` is the core radius and ``\alpha = \ln\left(\frac{\Delta r_{\max}}{\Delta r_{\min}}\right)``.

The nodal positions are then defined as:

```math
r_i = R_E + (R_c - R_E)
\left(
\frac{
\exp\left( \alpha \frac{N - i}{N - 1} \right) - 1
}{
\exp(\alpha) - 1
}
\right), \quad i = 1, \ldots, N
```

Rearranging the finite-difference equation yields a linear relation between adjacent nodes:

```math
\mathbf{C}_n \mathbf{y}_n + \mathbf{D}_{n+1} \mathbf{y}_{n+1} = \mathbf{0}
```

where the matrices ``\mathbf{C}_n`` and ``\mathbf{D}_{n+1}`` are defined as:

```math
\mathbf{C}_n = \mathbf{I} + \frac{\Delta r}{2} \mathbf{A}_n
```

```math
\mathbf{D}_{n+1} = -\mathbf{I} + \frac{\Delta r}{2} \mathbf{A}_{n+1}
```

The columns of ``\mathbf{I}_C`` represent the three elementary solutions of an isotropic, homogeneous elastic–gravitational sphere. The vector ``\mathbf{C}`` contains the integration constants:

```math
\mathbf{C} = (C_1, C_2, C_3)^T
```

These constants are determined from the boundary conditions at the planetary surface, applied to the stress components of the spheroidal solution. The general solution can be written as a linear combination of the elementary solutions ``\mathbf{y}^{(i)}``, which are regular at the centre:

```math
\mathbf{y}_l(r_C^+) = C_1 \mathbf{y}^{(1)} + C_2 \mathbf{y}^{(2)} + C_3 \mathbf{y}^{(3)}
```

Eliminating the coefficients ``C_i`` leads to:

```math
\begin{pmatrix}
U \\
V \\
\phi
\end{pmatrix}
-
\begin{pmatrix}
U^{(1)} & U^{(2)} & U^{(3)} \\
V^{(1)} & V^{(2)} & V^{(3)} \\
\phi^{(1)} & \phi^{(2)} & \phi^{(3)}
\end{pmatrix}
\begin{pmatrix}
X^{(1)} & X^{(2)} & X^{(3)} \\
Y^{(1)} & Y^{(2)} & Y^{(3)} \\
\psi^{(1)} & \psi^{(2)} & \psi^{(3)}
\end{pmatrix}^{-1}
\begin{pmatrix}
X \\
Y \\
\psi
\end{pmatrix}
= 0
```

Here, ``U^{(i)}, V^{(i)}, \phi^{(i)}, X^{(i)}, Y^{(i)}, \psi^{(i)}`` are the components of the (i)-th elementary solution.

This can be written compactly as a boundary condition:

```math
\mathbf{B} \, \mathbf{y} = 0
```

where ``\mathbf{B}`` is a ``3 \times 6`` matrix:

```math
\mathbf{B} =
\begin{pmatrix}
1 & 0 & b_{11} & b_{12} & 0 & b_{13} \\
0 & 1 & b_{21} & b_{22} & 0 & b_{23} \\
0 & 0 & b_{31} & b_{32} & 1 & b_{33}
\end{pmatrix}
```

The coefficients ``b_{ij}`` are defined as:

```math
\begin{pmatrix}
b_{11} & b_{12} & b_{13} \\
b_{21} & b_{22} & b_{23} \\
b_{31} & b_{32} & b_{33}
\end{pmatrix}
=
-
\begin{pmatrix}
U^{(1)} & U^{(2)} & U^{(3)} \\
V^{(1)} & V^{(2)} & V^{(3)} \\
\phi^{(1)} & \phi^{(2)} & \phi^{(3)}
\end{pmatrix}
\begin{pmatrix}
X^{(1)} & X^{(2)} & X^{(3)} \\
Y^{(1)} & Y^{(2)} & Y^{(3)} \\
\psi^{(1)} & \psi^{(2)} & \psi^{(3)}
\end{pmatrix}^{-1}
```

Notably, the second and third columns are swapped compared to Kobayashi (2007) for consistency with ``\mathbf{A}``. The difference equation together with the boundary conditions (given in the next sections) define the complete system of oscillations:

```math 
\begin{pmatrix}
B_1 \\
C_1 & D_2 \\
& C_2 & D_3 \\
& & \ddots & \ddots \\
& & & C_{N-1} & D_N \\
& & & & B_N
\end{pmatrix}
\begin{pmatrix}
y_1 \\
y_2 \\
y_3 \\
\vdots \\
y_{N-1} \\
y_N
\end{pmatrix}
=
\begin{pmatrix}
0 \\
0 \\
0 \\
\vdots \\
0 \\
b
\end{pmatrix}
```

Here, ``B_1`` and ``B_N`` represent the inner and outer boundary conditions, respectively, while ``C_n`` and ``D_n`` correspond to the discretised equations of motion, which depend on the forcing frequency ``\omega`` and the interior structure of the mantle. All unspecified entries in the global matrix are zero.

The structure of the system can be exploited to reduce computational cost. Following a Henyey-type relaxation method (e.g. Unno et al. 1989), the system can be rewritten in block form:

```math 
\begin{pmatrix}
S_1 & Q_1 \\
P_2 & S_2 & Q_2 \\
& \ddots & \ddots & \ddots \\
& & P_{N-1} & S_{N-1} & Q_{N-1} \\
& & & P_N & S_N
\end{pmatrix}
\begin{pmatrix}
y_1 \\
y_2 \\
\vdots \\
y_{N-1} \\
y_N
\end{pmatrix}
=
\begin{pmatrix}
0 \\
0 \\
\vdots \\
0 \\
b
\end{pmatrix}
```

In this formulation, ``S_1`` contains contributions from ``B_1`` and the upper block of ``C_1``. The upper and lower parts of ``Q_1`` consist of a zero block and the upper part of ``D_2``, respectively, while ``P_2`` contains the lower part of ``C_1`` combined with a zero block.

The submatrices are defined as:

```math 
P_n \equiv
\begin{bmatrix}
C^{l}_{n-1} \\
0
\end{bmatrix},
\quad n = 2, \ldots, N
```

```math
S_n \equiv
\begin{bmatrix}
D^{l}_n \\
C^{u}_n
\end{bmatrix},
\quad n = 2, \ldots, N-1
```

```math 
Q_n \equiv
\begin{bmatrix}
0 \\
D^{u}_{n+1}
\end{bmatrix},
\quad n = 1, \ldots, N-1
```

Here, superscripts ``u`` and ``l`` denote the upper and lower blocks of the corresponding matrices. The matrix ``S_N`` consists of the lower part of ``D_N`` combined with ``B_N``. Each of the block matrices ``P_n``, ``S_n``, and ``Q_n`` has size ``6 \times 6``. Both formulations are equivalent, and the full system corresponds to a banded matrix of size ``6N \times 6N`` for spheroidal oscillations in a solid sphere. Because coupling is restricted to nearest neighbours, the system can be solved efficiently using a recursive scheme.

The first block row reads:

```math
S_1 y_1 + Q_1 y_2 = 0
```

which gives:

```math 
y_1 = -S_1^{-1} Q_1 y_2
```

Substituting into the next block row yields:

```math 
\left(-P_2 S_1^{-1} Q_1 + S_2 \right) y_2 + Q_2 y_3 = 0
```

so that:

```math 
y_2 = -\left(-P_2 S_1^{-1} Q_1 + S_2 \right)^{-1} Q_2 y_3
```

Proceeding recursively, we obtain:

```math 
y_n = R_n y_{n+1}
```

with

```math 
R_n = -X_n^{-1} Q_n
```

and

```math 
X_n = P_n R_{n-1} + S_n
```

The initial condition is:

```math 
R_1 = -S_1^{-1} Q_1
```

Finally, the last block equation becomes (see below for the surface boundary conditions):

```math
X_N y_N = b =
\begin{pmatrix}
0 \\
0 \\
\frac{2n+1}{R}
\end{pmatrix}
```

with

```math 
X_N = P_N R_{N-1} + S_N
```

This recursive formulation reduces the computational cost from solving a full ``6N \times 6N`` system to a linear sequence of ``N`` matrix operations.

Two different propagation operators appear in the formulation. From the differential equations, one obtains the local transfer relation:

```math 
y_n = -C_n^{-1} D_{n+1} y_{n+1} \equiv \tilde{R}_n y_{n+1}
```

Although this resembles the recursive relation above, the operators are fundamentally different. The matrix ``\tilde{R}_n`` is constructed directly from the differential operator. Since both ``C_n`` and ``D_n`` are invertible, ``\tilde{R}_n`` is also invertible, allowing propagation in both directions along the grid.

In contrast, the relaxation operator is:

```math
y_n = R_n y_{n+1}, \qquad R_n = -X_n^{-1} Q_n
```

Here, ``Q_n`` is singular by construction, so ``R_n`` is not invertible. This enforces a one-way recursion from ``n+1`` to ``n``. Rather than being a limitation, this structure ensures numerical stability by enforcing simultaneous compatibility with both the differential equations and boundary conditions.


---

### Surface Boundary Conditions

Finally, to find the specific solution for the case where a planet is orbiting around a star-like object, we can specify the surface (``r=a^-``) boundary condition. In particular we need to constrain the 3rd, 4th, and 6th ``y``-values. The integration constants will follow from

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

To solve this system we thus need only provide ``\pmb{P}_1\,\pmb{y}(a^-)``. Additionally for the relaxation method we can impose a semi-free-surface outer boundary. The boundary conditions can be written in the same form as the inner system with

```math
\mathbf{B} =
\begin{pmatrix}
0 & 0 & 1 & 0 & 0 & 0 \\
0 & 0 & 0 & 1 & 0 & 0 \\
0 & 0 & 0 & 0 & (n+1)/R & 1
\end{pmatrix}
```

This corresponds to free-surface conditions applied to ``y_1``, ``y_2``, and (semi) ``y_5``, while the remaining components are constrained accordingly. In the general case, let ``P`` denote the external atmospheric pressure, ``\zeta`` the surface mass load, ``U`` an external potential, and ``\tau`` the tangential traction component. Then the boundary conditions for a spherical harmonic degree ``n`` can be written as:

```math
\begin{cases}
y_{3n}(a) = -g_e \zeta_n - P_n \\
y_{4n}(a) = \tau_n \\
y_{6n}(a) + \frac{n+1}{a} y_{5n}(a)
= \frac{2n+1}{a} U_n + 4\pi G \zeta_n
\end{cases}
```

where ``g_e`` is the surface gravity magnitude.

Farrell (1972) and Longman (1962) showed that the surface mass load can alternatively be expressed as an equivalent external potential ``U'``, such that:

```math
\zeta_n = \frac{2n+1}{4\pi G a} U'_n
```

We distinguish between Tidal Love numbers (response to an external potential perturbation ``U``) and Load Love numbers (response to a surface mass load ``\zeta`` (or equivalently ``U'``)). These can be obtained by solving the system with the following choices:

* TLN: ``(U, U', \tau, P) = (1, 0, 0, 0)``
* LLN: ``(U, U', \tau, P) = (0, 1, 0, 0)``

When the outer boundary is driven by tides from an external perturber, the final step in the relaxation method system can be written as:

```math
\mathbf{B} \mathbf{y} = \mathbf{b}
```

with

```math
\mathbf{b} =
\begin{pmatrix}
0 \\
0 \\
\frac{2n+1}{R}
\end{pmatrix}
```

This expression follows from a simplified form of the general boundary conditions, assuming the outer radius ``a = R`` (planetary radius).

---

