```@meta
CollapsedDocStrings = true
```

### Tutorials (1)

# Loading data

The `Obliqua` package comes with with several tidal heating modules. Each modules can be excessed through one cohensive function. Depending on which module you want to use `Obliqua` requires different input parameters. In order to get started, several data files are included, these can be used in combination with the different functions. First we shall show how to load these data files.

The data files are stored in the `/path/to/Obliqua/res` folder, and are of type `JSON`. They have following structure

```json
{
  "omega": "Float",
  "axial": "Float",
  "ecc": "Float",
  "sma": "Float",
  "S_mass": "Float",
  "density": "[Array]",
  "radius": "[Array]",
  "visc": "[Array]",
  "shear": "[Array]",
  "bulk": "[Array]",
  "bulkd": "[Array]",
  "phi": "[Array]",
  "perm": "[Array]",
  "cfg": "Dict"
}
```

Depending on which module is being used the following parameters need to be provided.

##### Table of inputs

| Input   | solid1d          | solid1d\_mush    | fluid1d  | Description                      | Symbol |
|:--------|:----------------:|:----------------:|:--------:|:---------------------------------|-------------:|
| omega   | ✔️              | ✔️               | ✔️       | Orbital Frequency                | ``\omega``   |
| axial   | ✔️              | ✔️               | ✔️       | Axial Frequency                  | ``\Omega``   |
| ecc     | ✔️              | ✔️               | ✔️       | Eccentricity                     | ``\epsilon`` |
| sma     | ✔️              | ✔️               | ✔️       | Semi major axis                  | ``a``        |
| S_mass  | ✔️              | ✔️               | ✔️       | Stellar mass                     | ``M_\star `` |
| density | ✔️              | ✔️               | ✔️       | Density profile                  | ``\rho``     |
| radius  | ✔️              | ✔️               | ✔️       | Radii                            | ``r``        |
| visc    | ✔️              | ✔️               | ✔️       | Viscosity profile                | ``\eta``     |
| shear   | ✔️              | ✔️               | ❌       | Shear profile                    | ``\mu``      |
| bulk    | ✔️              | ✔️               | ❌       | Solid Bulk Modulus               | ``\kappa_s`` |
| bulkd   | ❌              | ✔️               | ❌       | Drained Bulk Modulus             | ``\kappa_d`` |
| phi     | ❌              | ✔️               | ❌       | Melt Fraction profile (porosity) | ``\phi``     |
| perm    | ❌              | ✔️               | ❌       | Permeability profile             | ``k``        |
| cfg     | ✔️              | ✔️               | ✔️       | Configuration                    |  |

It is important to note that the `radius` array contains the radial values at the boundaries of the spherical shells that make up the planetary mantle, whilst the `density`, `visc`, `shear`, `bulk`, `bulk_d`, `phi`, and `perm` arrays contain the mean of these in the spherical shells. As such it follows that there are ``N+1`` values in the `radius` array, and ``N`` values in the `density`, `visc`, `shear`, `bulk`, `bulk_d`, `phi`, and `perm` arrays, where ``N`` is the number of spherical shells. 

For now, the user is provided with four different functions, they are given below. The most simple way to use these functions is as follows. First let's test if the provided data is compatible.

```julia
using Obliqua

# location of data files
RES_DIR = "/path/to/Obliqua/res"

# test data validity using included data file
ok = load.load_interior("$RES_DIR/interior_data/test_mantle.json", true)

if ok
    print("Pass")
else
    print("Fail")
end

```

Next, let's load the data from the data file.

```julia
# use the relevant load function
omega, ecc, rho, radius, visc, shear, bulk, ncalc = load.load_interior("$RES_DIR/interior_data/test_mantle.json", false)

```

Obliqua still expects the user to provide all quantities, in this case one can simply pass an array of zeros (e.g. for the `phi` array, since this is not used in the `solid1d` module.) Recommnende is to use 

```julia
omega, axial, ecc, sma, S_mass, rho, radius, visc, shear, bulk, phi, ncalc =
        load.load_interior_mush_full("$RES_DIR/interior_data/test_mantle_mush_full_test.json", false)
```

to get started, as it provides all the necessary parameters for all modules. The remaining parameters can be obtained using

```julia
cfg       = Obliqua.open_config("$RES_DIR/config/all_options.toml")

perm      = Obliqua.interior.get_permeability(phi, cfg)
perm, phi = Obliqua.interior.limit_porosity(perm, phi, cfg)
bulkd     = Obliqua.interior.get_drained_bulk(bulk, phi, cfg)
```

where `cfg` is the configuration file that contains all the necessary parameters for the tidal model. The `perm` and `bulkd` arrays are then calculated based on the `phi` and `bulk` arrays and the configuration file.

---
---

```@autodocs
Modules = [Obliqua.load]
Order   = [:function, :type]     
```

---
---