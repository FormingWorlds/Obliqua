
### Reference (1)

# Chapter 1 - Rheology

The frequency response of the mantle is computed for a homogeneous viscoelastic mantle. `fwlLove` provides two mantle rheologies, namely the Maxwell- and Andrade rheologies. The latter has been shown to be in good agreement with observational data. The rheologies are incoorporated through their definition of the complex shear modulus. For the Andrade rheology, the complex shear modulus is 

```math
\tilde{\mu}(\omega) =
\frac{\mu}{
1 + (\mathrm{i}\omega\eta/\mu)^{-\alpha}\Gamma(1+\alpha)
  + (\mathrm{i}\omega\eta/\mu)^{-1}}, \qquad \alpha = 0.3,
```

where ``\mu`` is the elastic shear modulus, ``\tau_A`` the Andrade timescale, ``\tau_M = \eta/\mu`` the Maxwell time, ``\alpha`` the Andrade exponent, and ``\Gamma`` the Gamma function.

For the Maxwell limit ``\alpha = 0`` and 

```math
\tilde{\mu}(\omega) =
\frac{\mu}{1 + (\mathrm{i}\omega\eta/\mu)^{-1}}.
```

---