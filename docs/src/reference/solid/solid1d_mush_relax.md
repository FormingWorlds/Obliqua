
### Reference (3)

# Solid-Phase - solid1d-mush-relax

For the "solid1d-mush-relax" model we will need to use many of the
concepts we saw for the other models. The basic structure is identical
to that of "solid1d-relax", but instead of just one motion matrix we
will now be working with two motion matrices (i.e. the two
$\pmb{A}_n(r)$ listed earlier in the document). In order to grasp what
is going on in the system let us draw schematically the structure of the
solver.

![Schematic of "solid1d-mush-relax"
model.](images/solid1d_mush_relax.png)

Note that for symmetry reasons I have reordered the $y$-functions. This
way when we embed the $6\times6$ system in the $8\times8$ system we
split the lower and upper parts and add one new porous $y$-function to
both sets. If we had not done this, then the new grouping would contain
4 old and 2 old + 2 new $y$-functions, instead of 3 old + 1 new and 3
old + 1 new. Within the code we create a transfer matrix array $R$
(length $N$) with size $8 \times 8$, we subsequently expose the
$6\times 6$ system as a "view". The specific ording as seen by the two
systems is 

$$\begin{aligned}
    \pmb{y}_{n,m} &= (U_{n,m} , V_{n,m} , \Phi_{n,m} , X_{n,m} , Y_{n,m}, \Psi_{n,m})^T & (6\times6) \\
    \pmb{y}_{n,m} &= (U_{n,m} , V_{n,m} , \Phi_{n,m} , P_{n,m},  X_{n,m} , Y_{n,m}, \Psi_{n,m}, R_{n,m})^T & (8\times8)
\end{aligned}$$ 

With this clarified, let us have a look at the diagram.

From the left, we start with the core solution boundary

$$
B =
\begin{pmatrix}
1 & 0 & 0 & b_{11} & b_{12} & b_{13} \\
0 & 1 & 0 & b_{21} & b_{22} & b_{23} \\
0 & 0 & 1 & b_{31} & b_{32} & b_{33}
\end{pmatrix},$$ 

We than impose continuity between $B$ and our modal
solution vector, this is denoted by the dotted vertical line.

$$\begin{pmatrix}
    1 & & & & &    \\
      & 1 & & &   \\
      & & 1 & &   \\
      & & & 1 &   \\
      & & & & 1   \\
      & & & & & 1  
    \end{pmatrix}_{(6\times6)}
    \pmb{y}_{n,m}(r_1^-) = 
    \begin{pmatrix}
    1 & & & & &    \\
      & 1 & & &   \\
      & & 1 & &   \\
      & & & 1 &   \\
      & & & & 1   \\
      & & & & & 1  
    \end{pmatrix}_{(6\times6)}
    \pmb{y}_{n,m}(r_C^+)$$ 
    
We then start constructing $R_n$ throughout
the solid layers 

$$
    \pmb{C}_n \pmb{y}_n + \pmb{D}_{n+1} \pmb{y}_{n+1} = \pmb{0}$$ 
    
with
$\pmb{A}_n$ the $6 \times 6$ motion matrix, exactly as in
[Solid-Phase - solid1d-relax](@ref). At some point we may encounter a
porous layer, where we must couple the $6\times6$ elastic system to the
$8\times8$ poro-elastic one that carries the two extra porous
$y$-functions (pore pressure $P_{n,m}$ and radial melt flux $W_{n,m} = \phi R_{n,m}$ relative to the matrix).
The melt flux is continuous wherever the pores connect, so a porosity change between two porous nodes needs no
interface condition.

### Solid-to-mush and mush-to-solid interfaces

At a transition, the response matrices are still built from the same
$\pmb{C}_n, \pmb{D}_{n+1}$ recursion, but now in the padded $8\times8$
space so that both sides of the interface can be expressed in a common
$8$-component vector. Two asymmetric adjustments are made to $\pmb C_n$
or $\pmb D_{n+1}$ depending on the direction of the transition, since one
side of the interface has no porous degrees of freedom to match against:

* **Solid $\to$ mush** (entering a porous layer from below): the incoming
  $3\times6$ "stored" lower half-block from the solid side is scattered
  into the $8$-slot ordering at the six non-porous column positions
  (i.e. everywhere except the pore-pressure and Darcy-flux slots), and
  the pore-pressure row of $\pmb C_n$ is zeroed, since the solid side has
  no pore pressure to enforce continuity on.
* **Mush $\to$ solid** (leaving a porous layer from below): both the
  pore-pressure and Darcy-flux rows of $\pmb D_{n+1}$ are zeroed, since
  the solid side above has neither degree of freedom.

At a **solid $\to$ mush** transition a small regularizing term $\pmb K_n$ (a single $1$ on the melt-flux
diagonal entry) is added to $\pmb X_n = \pmb P_n \pmb R_{n-1} + \pmb S_n$, $\pmb R_n = -\pmb X_n^{-1} \pmb Q_n$ is
obtained with a direct solve (the pseudo-inverse only if $\pmb X_n$ is singular), and the melt-flux row of $\pmb R_n$
is set to zero, so that $W = 0$ at the last solid node. The pivot only turns the (identically cancelled)
pore-pressure row of the interface cell into a slack; through the melt-flux row of that cell, $W$ then vanishes at
the cell face, which is the base of the porous layer (the pore pressure is free).

At a **mush $\to$ solid** transition the pore pressure jumps from its value on the porous side,
$(\pmb C_n \pmb y_n)_P$, to zero. With this self-referencing jump the pore-pressure row of the
difference equation across the interface vanishes identically, so the pore pressure is free. The
vacated row takes the melt-flux row of the same difference equation, $(\pmb C_n \pmb y_n)_W = 0$,
i.e. zero radial melt flux at the cell face (this row has no slot in the $6\times6$ solid block
row above), and $\pmb R_n = -\pmb X_n^{-1}\pmb Q_n$ is obtained with a direct solve. The $X$, $Y$, $\Psi$
rows of the interface cell keep their coupling to $P_n$ and $R_n$: that part of $\pmb C^l_n \pmb y_n$ is
returned as `Pn_extra` and added to $\pmb X_{n+1}$ by the next block row. With both choices the
interface cell is treated consistently and the scheme is second-order accurate across a buried porous
layer on a grid whose interfaces fall on cell faces. Imposing $R_n = 0$ at the last porous node instead is
first order (0.5 to 3 % error in $\mathrm{Im}\,k_2$ at $\Delta r = 1$ km in the test cases), and adding the
regularizing $\pmb K_n$ at this interface imposes $R_n = -(\pmb C_n \pmb y_n)_P$, a drained-type
condition that lets melt cross the interface (30 to 130 times too much Darcy dissipation in the test cases).

The Darcy dissipation is concentrated in a poroelastic boundary layer of thickness
$\delta = \sqrt{2 k M/(\eta_l \omega)}$, with $1/M = S + \alpha^2/(K_d + 4\mu/3)$, at every impermeable
boundary. On grids coarser than $\delta$, $\mathrm{Im}\,k_2$ remains close to its converged value, but the
heating of the unresolved boundary layer does not, and the energy budget fails (37 % at
$\Delta r = 20$ km for $\delta \approx 0.2$ km in the test case). See `test/test_mush_interfaces.jl`.

Once inside the porous layer, we continue constructing $R_n$ using the
full $8\times8$ motion matrix,

$$
    \pmb{C}_n \pmb{y}_n + \pmb{D}_{n+1} \pmb{y}_{n+1} = \pmb{0}, \qquad
    \pmb{A}_n \text{ the } 8\times8 \text{ motion matrix},$$

with the same forward/backward Henyey sweep used for the solid segments
above, until the mush layer ends (another interface of the opposite
kind) or the surface is reached.

A simple practical workaround is to assign a porosity slightly above the
percolation threshold (`porosity_thresh`) everywhere a porous layer is
present, so that the solver is built as a single, interface-free
$8\times8$ propagator instead of exercising this coupling logic at all.

---

### Function Documentation

```@docs
Obliqua.run_solid1d_mush_relax
```


### Heating profile

The cycle-averaged dissipation density is evaluated on the solver grid and integrated over the sphere on the angular grid
of `constants.res`. It is then volume-averaged onto the interior-model layers, which conserves the integrated heating.
There are three channels (Hay et al. 2025, Eqs. 38-39):

- **Shear:** $\omega\,\mathrm{Im}(\tilde\mu)\left[\sum_i |\epsilon_{ii}|^2 + 2\sum_{i<j}|\epsilon_{ij}|^2 - \tfrac13|\mathrm{tr}\,\epsilon|^2\right]$.
- **Darcy:** $\tfrac12\,\omega^2(\eta_l/k)\,|\boldsymbol w|^2$, with $\boldsymbol w = \phi\,\boldsymbol u_{rel}$ the melt flux relative to the matrix.
- **Compaction (volumetric):** with the solver's constitutive law (total mean stress $\tilde\kappa_d\,e - \tilde\alpha p$, fluid content $\tilde\alpha e + \tilde S p$, $e = \mathrm{tr}\,\epsilon$),

$$
\dot E_C = \frac{\omega}{2}\left[\mathrm{Im}(\tilde\kappa_d)\,|e|^2 - \mathrm{Im}(\tilde S)\,|p|^2 - 2\,\mathrm{Im}(\tilde\alpha)\,\mathrm{Re}(p\,e^*)\right]
         = \frac{\omega}{2}\,\mathrm{Im}(\tilde\kappa_d)\left|e + \frac{p}{\kappa_s}\right|^2,
$$

  where the second form uses $\tilde\alpha = 1 - \tilde\kappa_d/\kappa_s$ and $\tilde S = \phi/\kappa_l + (\tilde\alpha - \phi)/\kappa_s$ with real $\kappa_s, \kappa_l$.
  This follows from the work rate $\sigma_{ij}\dot\epsilon_{ij} + p\,\dot\zeta$ (the solid term of Hay et al. Eq. 38).
  Hay et al. Eq. 39b, $\tfrac{\omega}{2}\mathrm{Im}(\tilde\kappa_d)\,(|e|^2 + |p/\kappa_s|^2)$, lacks the cross term $2\,\mathrm{Re}(p\,e^*)/\kappa_s$.
  It is exact only for incompressible grains ($\kappa_s\to\infty$). Without the cross term, a porous layer with a lossy
  drained modulus overestimates the heating by ~9 % against the Love-number dissipation $-(2n+1)R\omega\,\mathrm{Im}(k_2)/(8\pi G)$.
  With it, the balance closes to < 0.1 % (`test/test_profile_coupling.jl`).
