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
  "phi": "[Array]",
  "ncalc": "Int"
}
```

Depending on which module is being used the following parameters need to be provided.

##### Table of inputs

| Input   | solid-phase      | solid-phase + mush interface | liquid-phase | Description | Symbol |
|:--------|:----------------:|:----------------:|:--------:|:---------------------------------|-------------:|
| omega   | ✔️              | ✔️               | ✔️       | Orbital Frequency                | ``\omega``   |
| axial   | ❌              | ❌               | ✔️       | Axial Frequency                  | ``\Omega``   |
| ecc     | ✔️              | ✔️               | ✔️       | Eccentricity                     | ``\epsilon`` |
| sma     | ❌              | ❌               | ✔️       | Semi major axis                  | ``a``        |
| S_mass  | ❌              | ❌               | ✔️       | Stellar mass                     | ``M_\star `` |
| density | ✔️              | ✔️               | ✔️       | Density profile                  | ``\rho``     |
| radius  | ✔️              | ✔️               | ✔️       | Radii                            | ``r``        |
| visc    | ✔️              | ✔️               | ✔️       | Viscosity profile                | ``\eta``     |
| shear   | ✔️              | ✔️               | ❌       | Shear profile                    | ``\mu``      |
| bulk    | ✔️              | ✔️               | ❌       | Solid Bulk Modulus               | ``\kappa_s`` |
| phi     | ❌              | ✔️               | ❌       | Melt Fraction profile (porosity) | ``\phi``     |
| ncalc   | ✔️              | ✔️               | ❌       | Resolution                       | ``n_\text{calc}`` |

It is important to note that the `radius` array contains the radial values at the boundaries of the spherical shells that make up the planetary mantle, whilst the `density`, `visc`, `shear`, `bulk`, and `phi` arrays contain the mean of these in the spherical shells. As such it follows that there are ``N+1`` values in the `radius` array, and ``N`` values in the`density`, `visc`, `shear`, `bulk`, and `phi` arrays, where ``N`` is the number of spherical shells. Please do not confuse `ncalc` with ``N``, since the former expands the number of layers through interpolation, whereas the latter is the initial resolution of the data. 

For now the user is provided with four different functions, they are given below. The most simple way to use these functions is as follows. First let's test if the provided data is compatible.

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

---
---

```@autodocs
Modules = [Obliqua.load]
Order   = [:function, :type]     
```

---
---