```@meta
CollapsedDocStrings = true
```

### Tutorials (2)

# Running Model

Now that you are able to validate and load data files using the `Obliqua.load` module, we can start using the tidal models  `Obliqua.solid0d`, `Obliqua.solid1d`, `Obliqua.solid1d_mush`, `Obliqua.solid1d_relax`, `Obliqua.solid1d_mush_relax`, `Obliqua.fluid0d` and `Obliqua.fluid1d`. In principle, this is rather simple, since you only have to call one function. For example, let us proceed with the data file from [Loading data](@ref). In this case, following the [Table of inputs](@ref), we are interested in the solid-phase tides only, so we may use the `Obliqua.calc_solid_tides` function. 

```julia
using Obliqua

# location of data files
RES_DIR = "/path/to/Obliqua/res"

# use the relevant load function
omega, axial, ecc, sma, S_mass, rho, radius, visc, shear, bulk, phi, ncalc =
    Obliqua.load.load_interior_mush_full("$RES_DIR/interior_data/test_mantle_mush_full_test.json", false)

# call the desired model
power_prf, power_blk, σ_range, imag_k2 = Obliqua.run_tides(
    omega, axial, ecc, sma, S_mass, rho, radius, visc, shear, bulk, phi, cfg
)

```

---
---


```@autodocs
Modules = [Obliqua]
Order   = [:function, :type]   
```

---
---