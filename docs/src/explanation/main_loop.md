
### Main loop

Now that the user is able to interface with Obliqua, we will shift our
focus to the inner workings of the main function. Before starting the
main iteration, Obliqua checks if all required configuration parameters
are specified. Subsequently, it loads the configuration into the
variables used by the code. Some configuration parameters may be set to
"none", the function `nothing_if_none` converts these strings to nothing
type elements in Julia. The interior arrays are then converted to
Bigfloat, this helps stabilize the solver and allows it to cover even
large discontinuities.

### Segments

The interior is subdivided into so-called "segments". These segments are
identified based on their local viscosity, any segment below the
liquidus viscosity is liquid, and above the solidus is solid. In between
these the formalism employed is mush. Notably, the advised choice it
setting the solidus viscosity to $\sim 10^{5}$ Pa s, allowing for
solid-like viscous dissipation and compaction to characterize the
dissipation in the mush. Only at very low viscosities does one expect
fluid waves, as such the liquidus may be set to as low as $\sim 10^2$ Pa
s. To bridge the region between these models one can either further
extend the solid formalism, or use the "interp" model in the mush,
allowing for a more smooth transition in the heating profile.

### Core

The code assumes a uniform core. The effects on the local gravity are
calculated and the density for the fluid density contrast are set.

### $k$-range

The tidal potential can be expanded in an infinite sum of spherical
harmonic terms with some degree $n$, order $m$, and harmonic $k$
(integer values). The user has the option the provide these quantities
in the configuration. However, specifically the $k$-range may also be
determined by the code. The desired range increases with eccentricity.
The Hansen coefficients $X_k = X^{-(n+1), m}_k(e)$ are computed by FFT
over the mean anomaly, for any eccentricity below 1 (see
[Forcing Frequency](@ref)). In the adaptive spectrum, unless `s_min` and
`s_max` are set, the range $K = [k_\mathrm{lo}, k_\mathrm{hi}]$ is chosen
by trimming modes from both ends of that spectrum.

Each mode $k$ gets the weight

```math
w_k = (1 + |k|)\, X_k^2 .
```

The two terms measure what a mode adds to the tide:

- $X_k^2$ is proportional to $|U_{n,m,k}|^2$, the square of the tidal
  potential of mode $k$ (see [Tidal potentials](@ref)).
- $|k| X_k^2$ tracks the heating of mode $k$. The heating is
  $H_k \propto \sigma_k\,(-\mathrm{Im}[k_n(\sigma_k)])\,X_k^2$ with
  $\sigma_k = m\Omega - k\,n_\mathrm{orb}$. The Love number is not known
  before the solve, so the weight assumes a frequency-independent phase
  lag (constant $\mathrm{Im}[k_n]$). Then $H_k \propto |\sigma_k| X_k^2$,
  and in the tails of the spectrum ($|k|\,n_\mathrm{orb} \gg m\Omega$)
  $|\sigma_k| \propto |k|$. The spin rate is left out, so the range
  depends only on $e$, $n$ and $m$.

With $W = \sum_k w_k$ the total weight and $\varepsilon$ the config value
`hansen_tail_tol` (default $10^{-6}$), $K$ is the smallest range for which
the weight left out on each side is below $\varepsilon W$:

```math
\sum_{k < k_\mathrm{lo}} w_k < \varepsilon\, W ,
\qquad
\sum_{k > k_\mathrm{hi}} w_k < \varepsilon\, W .
```

In practice modes are removed one at a time from the low end, then from
the high end, for as long as the removed weight on that side stays below
$\varepsilon W$. At most $2\varepsilon$ of the weight is dropped in
total. For $m = 0$ the spectrum is symmetric in $k$; only $k \geq 0$ is
kept, up to the larger of $|k_\mathrm{lo}|$ and $k_\mathrm{hi}$, and the
mirrored modes are counted twice. A circular orbit gives the single mode
$k = m$.

The cut-off bounds the weight, not the heating itself: if
$\mathrm{Im}[k_n]$ peaks at a frequency in the tail (a resonance), the
heating left out can exceed $\varepsilon$ of the total. A smaller
`hansen_tail_tol` widens $K$.

For testing convenience Obliqua includes two modes: "full" and
"adaptive". In principle, one should use "full" to test the tidal
calculation over a large range of forcing frequencies, and use
"adaptive" when using Obliqua with a PROTEUS-like framework. The former
will determine the planet response to all relevant forcing frequencies
and assumes $k = 1$ for every forcing frequency, while the latter only
calculates the response at the exact forcing frequencies excited by the
specific orbital configuration using $K$.

### Forcing frequency

The forcing frequencies $\sigma$ are simply obtained for every order $m$
and harmonic $k$ as 

```math
\sigma = m \Omega - k \omega.
```

### Complex shear modulus

The frequency response of the mantle is computed assuming a homogeneous
viscoelastic mantle. Obliqua implements two rheological models: Maxwell
and Andrade. Unlike the Maxwell model, which captures only
elastic--viscous behavior through a single Maxwell timescale, the
Andrade rheology includes a transient anelastic component. This yields
an elastic response at high forcing frequencies, a viscous response at
low frequencies governed by the Maxwell timescale $\tau_M = \eta / \mu$,
and a smooth anelastic transition in between characterized by the
Andrade timescale $\tau_A$ and creep exponent $\alpha$. As a result, the
Andrade rheology reproduces the observed attenuation behavior of
planetary mantles more accurately, particularly at high tidal forcing
frequencies. The rheologies are incorporated via their definitions of
the complex shear modulus. For the Andrade rheology, the complex shear
modulus is

```math
\tilde{\mu}(\sigma) =
\frac{\mu}{
1 + (\mathrm{i}\sigma\eta/\mu)^{-\alpha}\Gamma(1+\alpha)
  + (\mathrm{i}\sigma\eta/\mu)^{-1}}, \qquad \alpha = 0.3,
```

where $\mu$ is the elastic shear modulus and $\Gamma$ the Gamma
function.

For the Maxwell limit $\alpha = 0$ and

```math 
\tilde{\mu}(\sigma) =
\frac{\mu}{1 + (\mathrm{i}\sigma\eta/\mu)^{-1}}.
```
