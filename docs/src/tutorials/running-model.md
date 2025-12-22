```@meta
CollapsedDocStrings = true
```

### Tutorials (2)

# Running Model

Now that you are able to validate and load data files using the `Obliqua.load` module, we can start using the tidal models `Obliqua.Love` and `Obliqua.Fluid`. In principle, this is rather simple, since you only have to call one function depending on your use case. For example, let us proceed with the data file from [Loading data](@ref). In this case, following the [Table of inputs](@ref), we are interested in the solid-phase tides, so we should use `Obliqua.Love.calc_lovepy_tides`.

```julia
using Obliqua

# location of data files
RES_DIR = "/path/to/Obliqua/res"

# use the relevant load function
omega, ecc, rho, radius, visc, shear, bulk, ncalc = Obliqua.load.load_interior("$RES_DIR/interior_data/test_mantle.json", false)

# call the desired model
power_prf, power_blk, imag_k2 = Obliqua.calc_lovepy_tides(omega, ecc, rho, radius, visc, shear, bulk; ncalc=ncalc, material="maxwell")

```

In general, we can use 

| Tides                     | Function                               |
|:--------------------------|:---------------------------------------|
| solid-phase               | `Obliqua.calc_lovepy_tides`       |
| solid-phase + mush layer  | `Obliqua.calc_lovepy_tides_mush`  |
| liquid-phase              | `Obliqua.calc_fluid_tides`       |
| solid+liquid-phase        | `Obliqua.run_tides`                    |

---
---


```@autodocs
Modules = [Obliqua]
Order   = [:function, :type]   
```

---
---