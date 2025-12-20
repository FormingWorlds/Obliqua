
### Reference (4)

# Liquid-Phase

When the melt fraction of a magma layer exceeds a critical value ``F_m \gtrsim F_{m,c}``, its response to tidal forcing becomes fluid-like rather than solid-like. In this regime, tidal deformation is governed by the linearized Laplace tidal equations (LTEs), which describe the dynamics of a thin, global fluid shell of uniform thickness ``H`` on a planet of radius ``R_p``.

Neglecting nonlinear advection and mean flows, momentum and mass conservation for the magma ocean read
```math
\partial_t \vec{u} + \sigma_R \vec{u} + \vec{f} \times \vec{u} + g \nabla \zeta = g \nabla \zeta_{\mathrm{eq}}, \qquad
\partial_t \zeta + \nabla \cdot (H \vec{u}) = 0,
```
where ``\vec{u}`` is the horizontal velocity, ``\zeta`` the tidally induced surface displacement, ``\zeta_{\mathrm{eq}} = U_T/g`` the equilibrium tide associated with the tidal potential ``U_T``, ``g`` the surface gravity, and ``\sigma_R`` a Rayleigh drag frequency parameterizing dissipation in the magma ocean.

We focus on the strongly damped, highly viscous limit (``\sigma_R \gg 2\Omega``), appropriate for magma oceans, in which dissipative processes dominate over rotational effects. In this high-Ekman-number regime, the Coriolis term is negligible and the tidal response becomes overdamped. Dissipation is encapsulated by ``\sigma_R``, which should be interpreted as an effective damping rate accounting for viscous resistance, boundary-layer friction, form drag, and porous-flow dissipation near the rheological transition.

Transforming the governing equations to the frequency domain and expanding tidal quantities in spherical harmonics yields a closed-form expression for the complex, frequency-dependent tidal Love number ``k_\ell(\sigma)``. 

```math
k_{\ell}(\sigma)
= -\frac{\varrho_\ell \bar{\sigma}_\ell^2}{\sigma\tilde{\sigma} - \bar{\sigma}_\ell^2},
```

where ``\varrho_\ell = 3/(2n+1)(\rho_f/\rho_b)`` is the degree-``\ell`` density ratio between the fluid and the bulk part of the planet. The characteristic frequency is given by

```math
\bar{\sigma}_\ell = \sqrt{\frac{\mu_\ell g H_\mathrm{magma}}{R^2}},
```

where ``H_\mathrm{magma}`` the ocean depth, ``R`` the radius of the planet, and with 

```math
\mu_\ell = \ell(\ell+1),
```

the eigenvalue of the nth degree spherical harmonic, and 

```math
\tilde{\sigma} = \sigma - i \sigma_R.
```

---