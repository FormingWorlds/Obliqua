
### Reference (2)

# Rheology

The frequency response of each layer enters through its complex moduli: the shear modulus ``\tilde\mu(\sigma)`` and the (drained) bulk modulus. Obliqua uses the time dependence ``e^{\mathrm{i}\sigma t}``, so ``\mathrm{Im}\,\tilde\mu > 0`` for ``\sigma > 0`` is dissipation. The rheology is chosen separately for shear (`material_mu`) and bulk (`material_k`); for the bulk modulus replace ``\mu`` with ``\kappa`` and ``\eta`` with the bulk viscosity ``\zeta``, with compaction viscosity

```math
    \zeta (\phi) \approx \frac{\eta (\phi)}{\phi}.
```

Below, ``\mu`` is the unrelaxed (elastic) modulus, ``\eta`` the viscosity, ``\tau_M = \eta/\mu`` the Maxwell time and ``s = \mathrm{i}\sigma\tau_M``.

### Models in series

Maxwell, Andrade, Burgers and Sundberg-Cooper are elements in series, so their compliances ``J = 1/\tilde\mu`` add. Each starts from the elastic compliance ``1/\mu`` and the Maxwell dashpot ``1/(\mathrm{i}\sigma\eta)``:

```math
\frac{\mu}{\tilde\mu(\sigma)} = 1 + \frac{1}{s}
  \;+\; \underbrace{\Gamma(1+\alpha)\, s^{-\alpha}}_{\text{Andrade}}
  \;+\; \underbrace{\frac{\Delta}{1 + \mathrm{i}\sigma\tau_K}}_{\text{Kelvin-Voigt}} .
```

| `material` | Andrade term | Kelvin-Voigt term |
|:--|:--:|:--:|
| `"maxwell"` | – | – |
| `"andrade"` | ✔️ | – |
| `"burgers"` | – | ✔️ |
| `"sundberg-cooper"` | ✔️ | ✔️ |

- **Maxwell**: elastic at high frequency (``\sigma\tau_M \gg 1``), viscous at low frequency, with a single relaxation time ``\tau_M``.
- **Andrade** adds transient creep with exponent ``0 < \alpha < 1`` (`alpha`, default 0.3) and the Andrade time taken equal to ``\tau_M``. Its attenuation falls off as ``\sigma^{-\alpha}`` above the Maxwell frequency, as observed for mantle rocks, instead of the ``\sigma^{-1}`` of Maxwell.
- **Burgers** adds a Kelvin-Voigt element (a spring ``\mu/\Delta`` in parallel with a dashpot ``\eta_K``) in series with the Maxwell body. `kelvin_delta` sets ``\Delta``, the compliance of the element relative to ``1/\mu`` (default 0.2), and `kelvin_eta_ratio` sets ``\eta_K/\eta`` (default 0.02), which gives the relaxation time ``\tau_K = \Delta\eta_K/\mu``. Between the two times (``\tau_K \ll 1/\sigma \ll \tau_M``) the modulus sits on the relaxed plateau ``\mu/(1+\Delta)``, and the element adds a Debye peak to the loss compliance at ``\sigma\tau_K = 1``, of height ``\Delta/(2\mu)``.
- **Sundberg-Cooper** combines both terms: Andrade background with a Kelvin-Voigt peak (Sundberg & Cooper 2010), using the same three parameters.

All four are fully relaxed at ``\sigma = 0`` (``\tilde\mu = 0``). For ``\Delta = 0`` Burgers reduces exactly to Maxwell and Sundberg-Cooper to Andrade. Note that ``\alpha = 0`` is not the Maxwell limit: it doubles the dashpot term instead; Maxwell is the model without the Andrade term.

### Models in parallel and limits

- **Voigt** (`"voigt"`, Kelvin-Voigt solid): spring and dashpot in parallel, ``\tilde\mu = \mu + \mathrm{i}\sigma\eta``. Elastic at ``\sigma = 0``, increasingly viscous at high frequency.
- **Newton** (`"newton"`, Newtonian fluid): ``\tilde\mu = \mathrm{i}\sigma\eta``; ``\mu`` is not used.
- **Elastic** (`"elastic"`, equivalent to `"none"`): ``\tilde\mu = \mu``, no dissipation. This is useful as a dissipation-free limit for testing and for isolating the elastic part of the tidal response.

Every model satisfies ``\tilde\mu(-\sigma) = \overline{\tilde\mu(\sigma)}``, so negative forcing frequencies give the same dissipation.

---

### Function Documentation

```@docs
Obliqua.complex_modulus
```