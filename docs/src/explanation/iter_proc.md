
### Iterative Process

Now that all the required quantities and settings have been determined
and checked, we can move to the actual tidal calculation. Obliqua loops
over the segments, at each step it updates the mean density of the layer
below the current for the density contrast (or density ratio)
calculation. Additionally, specifically for the fluid models we include
an efficiency term to account for the reduced dissipation between the
fluid mantle and the fluid core. Note that at the moment an intermediate
mush segment forms, the efficiency is set back to 1. Within each segment
it loops over all the forcing frequencies.

For each forcing frequency the code calls the corresponding tides model
(e.g. solid, fluid, mush). For now these models basically act as sort of
black boxes where some interior properties are provided and some forcing
frequency is given, the model then returns the expected planet-wide
deformation, a.k.a. the $n$th degree Lovenumber $k_n(\sigma)$; the
planet-wide loading Lovenumber $k'_n(\sigma)$; and the normalized
heating profile in the segment. All the details regarding these models
will be given in the corresponding sections below. The currently
available tidal models are "solid0d", "solid1d", "solid1d_relax",
"solid1d_mush", "solid1d_mush_relax"; "fluid0d", "fluid1d"; "interp",
"none". For details see the Reference documantation.

The "interp" model requires knowledge of heating at both interfaces, as
such an additional code block is included to update the heating in the
"interp" region during the tidal calculation in the next segment after
the "interp" region.
