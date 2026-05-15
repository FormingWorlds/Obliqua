
### Post-processing

The final steps involve collecting the Lovenumbers and heating, and
applying the corresponding scalings to the normalized quantities.

The Lovenumbers can be combined using the following relation

```math 
k_n(\sigma) = k_{n}^{(\mathrm{solid})}(\sigma) + \bigl[1 + k'_{n}(\sigma)\bigr]\,k_{n}^{(\mathrm{fluid})}(\sigma).
```

Note that this relation implicitly assumes a solid interior covered by a
fluid magma ocean. For our purposes this assumption is fine. In the
special case where there is a subsurface fluid-like magma ocean, we may
include a correction (see Farhat 2025) due to a lithosphere. However, if
the subsurface magma ocean occurs deep in the mantle, it may be better
modeled using the solid formalism. More appropriate would be to use of
the "solid1d_mush" models as it accounts also for pore pressure in
partial melts. I have not yet encountered a case with a fluid-like deep
subsurface magma ocean. Despite this, it should be noted that in this
extreme scenario the model will likely not be able to properly represent
dissipation in the fluid, since there is currently no way to properly
combine the Lovenumbers.

Next we shall construct the global heating profile and bulk heating. For
the latter we require the imaginary part of the global tidal Lovenumber

```math
-\Im[k_n(\sigma)]
```

If earlier we had define the spectrum "full" then at this point we shall
assume that $k=1$ and determine the corresponding Hansen coefficient
(note that we have to make this assumption here since there is no
meaningful choice for $k$ since we picked an arbitrary forcing frequency
range.) 

$$X_1^{-(n+1),m}(e)$$ 

and normalization

```math
A_{n,m,1} = (2 - \delta_{m,0}\delta_{1,0}) (1 - \delta_{m,0}\delta_{1<0}) \sqrt{\frac{2 {\times 2\pi}}{2n+1}\frac{(n-m)!}{(n+m)!}} P_n^m(0) X_1^{-(n+1),m}(e).
```

The tidal potential is

$$U_{n,m,1} = \frac{GM}{a} \left(\frac{R}{a}\right)^n A_{n,m,1}(e)$$ 

The prefactor is $$\text{prefactor} \, = \frac{(2n + 1)R}{8πG} \sigma$$ The
normalized heating profile is then simply

$$H(r, \sigma) = \tilde{H}(r, \sigma) \times |U_{n,m,1}|^2$$ 

where $\tilde{H}(r, \sigma)$ is the unormalized heating profile returned by
the black box tides models. The global heating follows from

$$H(\sigma) = \text{prefactor} \times (-\Im[k_n(\sigma)]) \times |U_{n,m,1}|^2$$

If instead we had defined the spectrum "adaptive" then at this point we
will loop over the $(n,m,k)$-pairs and determine all heating based on
the most general expressions without making any assumptions (actually,
we still assume co-planarity, hence no obliquity tides :o). For every
pair we have 

$$X_k^{-(n+1),m}(e)$$ 

and normalization

$$A_{n,m,k} = (2 - \delta_{m,0}\delta_{k,0}) (1 - \delta_{m,0}\delta_{k<0}) \sqrt{\frac{2 {\times 2\pi}}{2n+1}\frac{(n-m)!}{(n+m)!}} P_n^m(0) X_k^{-(n+1),m}(e).$$

The tidal potential is

$$U_{n,m,k} = \frac{GM}{a} \left(\frac{R}{a}\right)^n A_{n,m,k}(e)$$ 

The prefactor is 

$$\text{prefactor} \, = \frac{(2n + 1)R}{8πG} \sigma$$ 

The normalized heating profile is then simply

$$H(r, \sigma) = \tilde{H}(r, \sigma) \times |U_{n,m,k}|^2$$ 

where $\tilde{H}(r, \sigma)$ is the unormalized heating profile returned by
the black box tides models. The global heating follows from

$$H(\sigma) = \text{prefactor} \times (-\Im[k_n(\sigma)]) \times |U_{n,m,k}|^2$$

Finally, in both spectrum cases, the heating rates can be integrated
across the forcing frequency domain 

$$H(r) = \int H(r, \sigma) d\sigma$$

and 

$$H = \int H(\sigma) d\sigma$$ 

(note the integral is actually just a
sum, given the discrete finite set of forcing frequencies.) The
resulting heating profile and global heating rates can be (and should
be!) compared through integration of the radial heating profile. They
should match as long as the surface load Lovenumber
$k'_n(\sigma) \ll 1$. For coupling to orbital dynamics we also return
the tidal Lovenumbers $k_n(\sigma)$ and corresponding forcing
frequencies $\sigma$.
