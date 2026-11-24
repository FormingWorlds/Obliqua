
### Reference (3)

# Solid-Phase - fluid1d_relax

`fluid1d_relax` (`module_fluid = "fluid1d_relax"`) solves the whole interior as one non-rotating, radially resolved column. Solid, mush and fluid layers, in any order, are passed through a single Henyey relaxation from the core-mantle boundary (CMB) to the surface, so tilt, loading and self-attraction of a fluid layer follow from the column solution rather than from separate terms. Without rotation the spherical-harmonic degrees decouple: only the forcing degree $n$ is solved, and the Love number is read at the surface,

$$k_n = y_5(R) - 1 .$$

---

### Layer classes

| Class | States | Motion matrix |
|---|---|---|
| solid | 6 | viscoelastic $6\times6$ matrix with the complex moduli of the layer (as in [Solid-Phase - solid1d-relax](@ref)) |
| mush | 8 | poro-viscoelastic $8\times8$ matrix (as in [Solid-Phase - solid1d-mush-relax](@ref)) |
| fluid | 6 | $6\times6$ matrix with the layer's own complex shear modulus ($\mu_c \approx i\omega\eta$ for a melt), the undrained (Gassmann) bulk modulus, and Rayleigh drag |

The fluid layers are the layers connected to the surface that satisfy `fluid_criterion`:

* `:viscosity`: viscosity below `visc_fluid`.
* `:skeleton`: the drained frame modulus is zero (failed skeleton). The drag then also acts on the solid and mush layers below, down to where the loss modulus stops exceeding the storage modulus at the forcing frequency.

Rayleigh drag acts on the horizontal motion only. In the tangential momentum row,

$$\rho\,\omega^2 \;\rightarrow\; \rho\,\omega\,(\omega - i\sigma_R).$$

By default (`profile = "nodes"`) the interior layer values are interpolated between the layer centres to every node of the solver grid, so the material and the fluid interface vary continuously. With `profile = "layers"` the values are constant per input layer, or per sub-layer with `n_sub > 1`.

In a mush the porous state is the melt flux $W = \phi R$, the volume flux of melt relative to the matrix ($R$ the relative melt displacement). The fluid-content equation conserves it, so $W$ and the pore pressure $P$ are continuous wherever the pores connect, and a change of porosity needs no interface condition (with $R$ as state, every porosity change would need the condition $\phi^- R^- = \phi^+ R^+$).

---

### State vector

The solver orders the states so that the upper half holds displacement-like quantities and the lower half the tractions:

| Slot | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 |
|---|---|---|---|---|---|---|---|---|
| Component | $U$ | $V$ | $\Phi$ | $P$ | $X=\sigma_{rr}$ | $Y=\sigma_{r\theta}$ | $\Psi$ | $W$ |

Solid and fluid nodes carry the six slots $\{1,2,3,5,6,7\}$; mush nodes carry all eight. $P$ is the pore pressure and $W$ the radial melt flux relative to the matrix.

---

### Interfaces and core

Every face is a set of linear conditions $\pmb L\,\pmb y^- = \pmb G\,\pmb y^+$ on the face values $\pmb y^- = \pmb C_n \pmb y_n$ and $\pmb y^+ = -\pmb D_{n+1}\pmb y_{n+1}$, which gives the cell equation

$$(\pmb L \pmb C_n)\,\pmb y_n + (\pmb G \pmb D_{n+1})\,\pmb y_{n+1} = \pmb 0 .$$

| Face | Condition |
|---|---|
| solid/mush $\leftrightarrow$ fluid, liquid core | free slip: $Y = 0$ on the side(s) that carry shear, $V$ free |
| solid $\leftrightarrow$ mush | closed pores: $W = 0$ on the porous side, $P$ free |
| mush $\leftrightarrow$ fluid | open pores: $U_f = U_s + W$, $P + X_f = 0$ |
| all other faces | all states continuous |

`slip = :noslip` and `pores = :closed` switch the fluid faces to the alternatives. The core enters through the left null space of its solutions $\pmb I_C$ (`get_Ic`); a liquid core's $\pmb I_C$ carries the tangential slip as its own column (free slip). A liquid core is solved as inertial-liquid, and as the static liquid core where the inertial solution is not finite (forcing frequencies below about $10^{-11}$ rad s$^{-1}$, where the two agree). A zero-frequency mode is the equilibrium tide (`solid1d-equil-relax`).

---

### Grid and surface

Nodes sit at the cell centres of Obliqua's exponential grid (`dr_min`, `dr_max`; see [Solid-Phase - solid1d-relax](@ref)). The CMB rows act on the face value $\pmb y(r_C) = (\pmb I - h_1 \pmb A_1)\,\pmb y_1$ and the surface rows on $\pmb y(R) = (\pmb I + h_N \pmb A_N)\,\pmb y_N$. Gravity is evaluated at the nodes from the enclosed mass. With these choices $k_n$ converges at second order in $\Delta r$.

---

### Heating and acceptance checks

The heating of the same solution has four channels: shear, Biot bulk, Darcy, and drag. Drag heating is

$$\dot e_d = \tfrac12\,\rho\,\sigma_R\,\omega^2\,|\pmb u_h|^2 ,$$

and in the fluid the strain is taken from the displacements, since $y_4/\mu_f$ is ill conditioned for $\mu_f \to 0$. Each solve is checked against:

* passivity: $-\mathrm{sign}(\omega)\,\mathrm{Im}\,k_n \ge 0$;
* the energy theorem
  $$\frac{(2n+1)\,R\,\omega}{8\pi G}\,\bigl(-\mathrm{Im}\,k_n\bigr) = \int \dot e \, dV ,$$
  with a warning when the ratio departs from 1 by more than `closure_tol`;
* resolution of the poroelastic boundary layers in the mush cells.

The heating profile is volume-averaged from the solver cells back onto the input layers, which conserves the integrated heating. The column model needs no load Love number, so it returns $k_L = 0$.

---

### Function Documentation

```@docs
Obliqua.run_fluid1d_relax
```
