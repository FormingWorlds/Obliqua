
### Reference (5)

# Forcing Frequency

Given the fact that the tidal forcing magnitude decreases exponentially with harmonic order, we may limit the calculation to only the lowest harmonic degree ``n = 2``. Generally, in the limit of ``R_p \ll a``, it suffices to only consider the quadrupolar harmonic (``n = 2``). Nevertheless, Obliqua can also determine higher degree contributions. As of now, we are considering a coplanar geometry. Hence, terms with ``m = 1`` associated with obliquity tides vanish. The Fourier modes of the second harmonic have frequencies given by

```math
\sigma = m\Omega - k n_{\mathrm{orb}},
```

where ``\Omega`` is spin rate and ``n_{\mathrm{orb}}`` orbital mean motion, and for integer values of order ``-2 \leq m \leq 2`` and harmonic ``\infty \leq k \leq \infty``. In our formalism tides are occuring over a large time interval, a time step ``\Delta t``. As such, we must account for tidal excitations that occur over a wide range of frequencies. We calculate the imaginary part of the ``n``th harmonic degree (``k_n``) Love number (``\Im[k_{n}(\sigma)]``) for all relevant harmonnic frequencies for which the Hansen coefficient 

```math
X_{nmk}(e) =
\frac{1}{2\pi} \int_0^{2\pi}
\left(\frac{r}{a}\right)^n
e^{im\Omega - ikn_{\mathrm{orb}}}\,dn_{\mathrm{orb}} \geq 0.01,
```

(i.e. ~1% corrections). 


---