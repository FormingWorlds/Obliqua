
### Reference (4)

# Mush layer - interp

The `interp` model is a simplified approach to modeling dissipation within a planetary region where active tidal equations (like those in `solid1d`) are not explicitly solved, such as a thin mushy layer or a transition zone. Instead of solving the full poro-viscoelastic system, it approximates the heating profile and subsequent Love numbers using exponential decay from the layer boundaries.

### Methodology

The model assumes that tidal dissipation peaks at the interfaces (upper and lower) and decays exponentially into the interior of the segment.

#### 1. Heating Profile Generation
The heating profile $P(r)$ is constructed as a superposition of two exponential decay functions originating from the top ($r_{top}$) and bottom ($r_{bot}$) interfaces:

$$P(r) = P_t \exp\left( -\frac{|r - r_{top}|}{l_t} \right) + P_b \exp\left( -\frac{|r - r_{bot}|}{l_b} \right)$$

where:
*   **$P_t, P_b$**: The heating intensities at the top and bottom interfaces.
*   **$l_t, l_b$**: The characteristic decay lengths, defined as a fraction of the total segment thickness ($l = \text{width} \cdot \Delta R$).



#### 2. Inferring Love Numbers
To remain consistent with the energy dissipation within the global model, the imaginary part of the Tidal Love number $k_2^T$ is inferred directly from the generated power profile. For a shell with volume $V_s$, the contribution to the complex Love number is:

$$\text{Im}(k_{2,s}^T) = -\frac{5 R \omega}{8\pi G} (P_s V_s)$$

The total $k_2^T$ for the segment is the sum of these contributions across all radial shells. Note that in this interpolation mode, the load Love number $k_2^L$ is typically returned as zero (or ignored) as the model focuses purely on energy dissipation rather than structural loading response.

---

### Function Documentation

```@docs
Obliqua.run_interp
```

### Configuration Parameters
When using this model via the TOML configuration, the following parameters in `[orbit.obliqua.mushy]` are relevant:

*   `t_width`: Controls the decay length from the top interface ($l_t$).
*   `b_width`: Controls the decay length from the bottom interface ($l_b$).

This model is particularly useful for representing "mushy" regions where the physical properties are highly uncertain, but the dissipation is expected to be concentrated near the boundaries of solid or liquid layers.