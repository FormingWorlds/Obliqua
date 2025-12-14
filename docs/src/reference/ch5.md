
### Reference (5)

# Chapter 5 - Forcing Frequency

The tidal power dissipation follows *Farhat et al. (2025, Eq. 30)*. The imaginary part (``\Im[k_{22}(\sigma)]``) describes the phase lag and is interpolated across positive and negative frequencies.

Each Fourier mode has frequency

```math
\sigma = 2\Omega - k n_{\mathrm{orb}},
```

where ``\Omega`` is spin rate and ``n_{\mathrm{orb}}`` orbital mean motion. The mode amplitude of the external tidal potential is

```math
U_{22k} =
\frac{G M_\star}{a}
\left(\frac{R}{a}\right)^2
\sqrt{\frac{6\pi}{5}}\,X_{22k}(e),
```

where ``X_{22k}(e)`` is the Hansen coefficient:

```math
X_{22k}(e) =
\frac{1}{2\pi} \int_0^{2\pi}
\left(\frac{r}{a}\right)^2
e^{2\mathrm{i}v - \mathrm{i}kM}\,dM.
```

Tidal power per mode:

```math
P_{T,k} =
\frac{5 R \sigma}{8\pi G}
\,\Im\!\bigl[k_{22}(\sigma)\bigr]
\,|U_{22k}|^2.
```

Total tidal dissipation:

```math
P_{\mathrm{tidal}} = -\sum_k P_{T,k}.
```