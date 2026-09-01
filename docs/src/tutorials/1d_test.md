

# One-dimensional test case: The Earth

This is an advanced test case to demonstrate the use of Obliqua for a one-dimensional model of the Earth. We will create a basic data lookup table for the Earth's tidal ``k``-Love number response, and map out the radial distribution of thedissipative tidal response. 

---

### Prerequisite

- Ensure that you have the `Obliqua` package installed.
- Attempt the zero-dimensional test case first.

---

### Creating the data file

This time the interior data will be provided, so it can simply be loaded using the `Obliqua.load` module. The file we will use is `res/interior_data/channel_showcase_updated.json`. 

---

### The configuration file

The configuration file should be a TOML file with the following structure:

```toml
title = "Earth"
version = "1.0"

[params]
    [params.out]
        path        = "out"
        time        = 0.0
        logging     = "INFO"
        plot_fmt    = "png"

# Orbital system
[orbit]
    [orbit.obliqua]
        store_3D    = false
        enforce_ec  = true
        optimize_scales = false
        solid_shell = true

        min_frac    = 0.05
        visc_l      = 1e2
        visc_lus    = 5e4
        visc_s      = 1e22
        visc_sus    = 5e6
        n           = [2]
        m           = [0, 2]
        spectrum    = "full"
        N_sigma     = 100
        p_min       = -8
        p_max       = 4
        s_min       = "none"
        s_max       = "none"
        
        material_mu = "andrade"
        material_k  = "andrade"
        alpha       = 0.3


        module_solid = "solid1d-relax"
        module_mushy = "interp"
        module_fluid = "fluid1d"

        [orbit.obliqua.solid]
            ncalc       = 100
            dr_min      = 10
            dr_max      = 3000
            core        = "liquid"
            core_props  = "core"
            inertial_terms = true
            bulk_l      = 1e9
            dbulk_power = 0.5
            porosity_thresh = 3e-2

        [orbit.obliqua.fluid]
            sigma_R     = 2e-3
            sigma_R_inf = 1e-3
            sigma_R_prf = "dynamic_interp"
            H_R         = 1e5
            efficiency  = 1e0

        [orbit.obliqua.mushy]
            b_width     = 5e-1
            t_width     = 3e-2

[struct]
    core_density = 10738.33
    core_shear   = 0.0
    core_bulk    = 5e11

[interior_energetics]
    grain_size   = 0.1
```

Create it in the `res/config` directory of the Obliqua package and name it `earth_config.toml`. Key parameters to note are the `module_solid`, `module_mushy`, and `module_fluid` parameters, which specify the tidal model to use for the solid, mushy, and fluid layers, respectively. In this case, we are using the `solid1d_relax` model for the solid layer and the `fluid1d` model for the fluid layer. We are also using the `andrade` rheological model for both the shear and bulk moduli. For the full spectrum, we are using 100 probe forcing frequencies, which are logarithmically spaced between 10⁻⁸ and 10⁴ Hz. 


### Running the test case

With these files in place, we can now run the test case. The following code snippet demonstrates how to do this.

```julia
using Obliqua

# Load configuration
cfg = Obliqua.open_config("res/config/earth_config.toml")

# Load interior model
omega, axial, ecc, sma, S_mass, rho, radius, visc, shear, bulk, phi, ncalc =
    load.load_interior_mush_full("res/interior_data/channel_showcase_updated.json", false)

# Extract mush zone properties
perm      = Obliqua.interior.get_permeability(phi, cfg)
perm, phi = Obliqua.interior.limit_porosity(perm, phi, cfg)
bulkd     = Obliqua.interior.get_drained_bulk(bulk, phi, cfg)

# Run tidal calculation
power_prf, power_blk, nmk, σ_range, LNk = Obliqua.run_tides(
    omega, axial, ecc, sma, S_mass, rho, radius, visc, shear, bulk, bulkd, phi, perm, cfg
)
```

---

### Expected output

At this point, the function will return several outputs in the terminal:

```bash
[ Info: Loading interior from JSON file: res/interior_data/channel_showcase_updated.json
[ Info: Using configuration 'Earth'
[ Info: Using (n, m, k) = (2, 0, 1) for full spectrum.
[ Info: Smoothing complex modulus profiles to avoid sharp jumps in viscosity...
[ Info: Smoothing between indices 15 and 20 (r = 3.14531757e6 to 3.21782745e6 m)
[ Info: Writing results to /home/marijn/LovePy/fwlLove.jl/out/0_obliqua.nc
[ Info: Expected bulk heating: 2.4659699632690844e14
[ Info: Obtained bulk heating: 2.4659699632690838e14
```

First, Obliqua will load the interior data and the configuration file. Next it will decide on which modes to include in the calculation. In this case, we are using the full spectrum, which restricts to the same mode for all forcing frequencies. Next, the complex modulus profiles are smoothed to avoid sharp jumps in viscosity, which can cause numerical instabilities. Finally, the results are written to a NetCDF file in the `out` directory. The expected bulk heating is also printed to the terminal, which is the total tidal dissipation rate in Watts. However, this is not relevant for the full spectrum. The obtained bulk heating is also printed, which is the total radially integrated tidal dissipation rate from the heating profile in Watts, and should match the expected value for all one-dimensional models.

Let us now run the `examples/visualize_output.ipynb` notebook to visualize the results. The output should look like this:

![alt text](panels/1d.png)

Here we see the imaginary part of the tidal ``k``-Love number as a function of the forcing frequency. We see both the solid and fluid contributions to the tidal response. The solid contribution is dominant at low frequencies, while the fluid contribution is dominant at high frequencies. This is the expected behaviour for a partially molten Earth. In addition, we observe high frequency resonances in the tidal response.




