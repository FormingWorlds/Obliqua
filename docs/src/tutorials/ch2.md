```@meta
CollapsedDocStrings = true
```

### Tutorials (2)

# Chapter 2 - Running Model

Now that you are able to validate and load data files using the `fwlLove.load` module, we can start using the tidal models `fwlLove.Love` and `fwlLove.Fluid`. In principle, this is rather simple, since you only have to call one function depending on your use case. For example, let us proceed with the data file from [Chapter 1 - Loading data](@ref). In this case, following the [Table of inputs](@ref), we are interested in the solid-phase tides, so we should use `fwlLove.Love.calc_lovepy_tides`.

```julia

using fwlLove

# location of data files
RES_DIR = "/path/to/fwlLove/res"

# use the relevant load function
omega, ecc, rho, radius, visc, shear, bulk, ncalc = load.load_interior("$RES_DIR/interior_data/test_mantle.json", false)

# call the desired model
power_prf, power_blk, imag_k2 = Love.calc_lovepy_tides(omega, ecc, rho, radius, visc, shear, bulk; ncalc=ncalc, material="maxwell")

```

In general, we can use 

| Tides                     | Function                               |
|:--------------------------|:---------------------------------------|
| solid-phase               | `fwlLove.calc_lovepy_tides`       |
| solid-phase + mush layer  | `fwlLove.calc_lovepy_tides_mush`  |
| liquid-phase              | `fwlLove.calc_fluid_tides`       |
| solid+liquid-phase        | `fwlLove.run_tides`                    |

---
---


```@autodocs
Modules = [fwlLove]
Order   = [:function, :type]   
```

---
---