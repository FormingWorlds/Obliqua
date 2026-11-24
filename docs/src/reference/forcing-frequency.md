
### Reference (1)

# Forcing Frequency

Given the fact that the tidal forcing magnitude decreases exponentially with harmonic order, we may limit the calculation to only the lowest harmonic degree ``n = 2``. Generally, in the limit of ``R_p \ll a``, it suffices to only consider the quadrupolar harmonic (``n = 2``). Nevertheless, Obliqua can also determine higher degree contributions. As of now, we are considering a coplanar geometry. Hence, terms with ``m = 1`` associated with obliquity tides vanish. The Fourier modes of the second harmonic have frequencies given by

```math
\sigma = m\Omega - k n_{\mathrm{orb}},
```

where ``\Omega`` is spin rate and ``n_{\mathrm{orb}}`` orbital mean motion, and for integer values of order ``-2 \leq m \leq 2`` and harmonic ``-\infty \leq k \leq \infty``. In our formalism tides are occuring over a large time interval, a time step ``\Delta t``. As such, we must account for tidal excitations that occur over a wide range of frequencies. We calculate the imaginary part of the ``n``th harmonic degree (``k_n``) Love number (``\Im[k_{n}(\sigma)]``) for all relevant harmonic frequencies. The Hansen coefficients

```math
X^{-(n+1), m}_k(e) =
\frac{1}{2\pi} \int_0^{2\pi}
\left(\frac{r}{a}\right)^{-(n+1)}
e^{imv}\, e^{-ikM}\,dM ,
```

with ``v`` the true and ``M`` the mean anomaly, are computed by FFT over ``M`` (Kepler's equation solved by a bracketed Newton iteration, so any ``e < 1`` is supported). In the adaptive spectrum the set ``K = [k_\mathrm{lo}, k_\mathrm{hi}]`` of harmonics is chosen from this spectrum: with weights ``w_k = (1 + |k|) X_k^2``, ``K`` is the smallest range that leaves less than ``\varepsilon\sum_k w_k`` outside on each side, where ``\varepsilon`` is the config value `hansen_tail_tol` (default ``10^{-6}``). ``X_k^2`` measures the tidal potential of mode ``k`` and ``|k| X_k^2`` its heating under a frequency-independent phase lag. A configured `s_min` or `s_max` replaces the matching end. For ``m = 0`` the spectrum is symmetric in ``k`` and only ``k \geq 0`` is kept. The weight is explained under "k-range" on the [main loop page](../explanation/main_loop.md).

---

### Function Documentation

```@docs
Obliqua.Hansen.hansen_modes
Obliqua.Hansen.get_hansen
Obliqua.Hansen.hansen_spectrum
Obliqua.Hansen.kepler_safe
```