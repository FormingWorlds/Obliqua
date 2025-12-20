
### Reference (7)

# Tidal potentials

Different sources can provide different tidal potentials, first we need to complete the solid-phase by providng the relevant surface boundary condition. Here we list some examples.

### Loaded Surface

For an impermeable load,

```math
y_2 =  - \frac{(2\ell + 1) g(R)}{4 \pi G R}, \qquad y_4 = 0, \qquad y_6 = \frac{2\ell + 1}{R}, \qquad y_8 = 0
```

where ``g(R)`` is the gravitational acceleration at the surface (e.g., Saito, 1974).


### Lunar Tide

For a fully solid mantle, the boundary conditions at the surface are

```math
y_3 = y_4 = 0, \qquad y_6 = - \frac{(2\ell + 1) g_s}{4} \frac{M_M}{M_E} \left(\frac{R_E}{a}\right)^3.
```

(see *Takeuchi et al. 1972*).

The final condition is

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

After having solved and generated the ``k_2 (\sigma)`` spectrum, we can find the power input from tides.

### Tidal Dissipation

The mode amplitude of the external tidal potential is

```math
U_{22k} =
\frac{G M_\star}{a}
\left(\frac{R}{a}\right)^2
\sqrt{\frac{6\pi}{5}}\,X_{22k}(e),
```

The tidal power per mode is

```math
P_{T,k} =
\frac{5 R \sigma}{8\pi G}
\,\Im\!\bigl[k_{2}(\sigma)\bigr]
\,|U_{22k}|^2.
```

Such that the total tidal dissipation:

```math
P_{\mathrm{tidal}} = -\sum_k P_{T,k}.
```

---