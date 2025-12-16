
### Reference (5)

# Forcing Frequency

Given the fact that the tidal forcing magnitude decreases exponentially with harmonnic order, we limit the calculation to only the lowest harmonnic frequency ``n = 2``. Generally, in the limit of ``R_p = a``, it suffices to only consider the quadrupolar harmonic (``n = 2``). Moreover, as we are considering a coplanar geometry, terms with ``\ell = 1`` associated with obliquity tides vanish. The Fourier modes of the second harmonnic have frequencies given by

```math
\sigma = 2\Omega - k n_{\mathrm{orb}},
```

where ``\Omega`` is spin rate and ``n_{\mathrm{orb}}`` orbital mean motion, and for integer values of ``k``. In our formalism tides are occuring over a large time interval, a time step ``\Delta t``. As such, we must account for tidal excitations that occur over a wide range of frequencies. We calculate the imaginary part of the second harmonnic degree (``k_2``) Love number (``\Im[k_{2}(\sigma)]``) for a handful of frequencies. The resulting descrete profile can subsequently be interpolated across positive and negative frequencies to yield a continous profile. 

Note, since we are only interested in the the Fourier modes of the second harmonnic, we don't actually use the entire spectrum. We can determine which frequencies are relevant using the Hansen Coefficients. 

```math
X_{22k}(e) =
\frac{1}{2\pi} \int_0^{2\pi}
\left(\frac{r}{a}\right)^2
e^{2i v - ikM}\,dM,
```

where ``v`` is the true anomoly. 

---