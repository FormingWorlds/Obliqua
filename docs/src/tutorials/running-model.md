```@meta
CollapsedDocStrings = true
```

### Tutorials (3)

# Running Model

Now that you are able to validate and load data files using the `Obliqua.load` module, we can start using the tidal models  `Obliqua.solid0d`, `Obliqua.solid1d`, `Obliqua.solid1d_mush`, `Obliqua.solid1d_relax`, `Obliqua.solid1d_mush_relax`, `Obliqua.fluid0d` and `Obliqua.fluid1d`. In principle, this is rather simple, since you only have to call one function. For example, let us proceed with the data file from [Loading data](@ref). 

```julia
using Obliqua

# location of data files
RES_DIR = "/path/to/Obliqua/res"

# use the relevant load function
omega, axial, ecc, sma, S_mass, rho, radius, visc, shear, bulk, phi, ncalc =
    Obliqua.load.load_interior_mush_full("$RES_DIR/interior_data/test_mantle_mush_full_test.json", false)

# Extract mush zone properties
perm      = Obliqua.interior.get_permeability(phi, cfg)
perm, phi = Obliqua.interior.limit_porosity(perm, phi, cfg)
bulkd     = Obliqua.interior.get_drained_bulk(bulk, phi, cfg)

# call the desired model
power_prf, power_blk, nmk, σ_range, LNk = Obliqua.run_tides(
    omega, axial, ecc, sma, S_mass, rho, radius, visc, shear, bulk, bulkd, phi, perm, cfg
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