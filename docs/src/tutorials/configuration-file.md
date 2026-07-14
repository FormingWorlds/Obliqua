```@meta
CollapsedDocStrings = true
```

### Tutorials (2)

# Configuration file

You can read the configuration from a TOML file using the following command:

```julia
cfg = Obliqua.open_config("$RES_DIR/config/all_options.toml")
```

The configuration files follow the conventions used within PROTEUS. The default `all_options.toml` file contains all available parameters.

### Global Parameters
These define metadata and output behavior:
* **`title`**: Identifier for the simulation setup.
* **`version`**: Configuration file version for reproducibility.

---

### Execution Parameters (`[params]`)
This block controls output and logging:
* **`path`**: Directory where output files are stored.
* **`time`**: Current time of the simulation run in years (used for output file naming).
* **`logging`**: Logging level (e.g., `INFO`, `DEBUG`).
* **`plot_fmt`**: Output format for generated plots.

---

### Stellar Parameters (`[star]`)
Defines the host star:
* **`mass`**: Stellar mass in solar masses ($M_\odot$).

---

### Planetary Parameters (`[planet]`)
Defines the planet's physical properties:
* **`mass_tot`**: Total planetary mass in Earth masses ($M_\oplus$).

---

### Orbital Parameters (`[orbit]`)
Describes the planetary orbit and optional satellite:
* **`semimajoraxis`**: Orbital semi-major axis [AU].
* **`eccentricity`**: Orbital eccentricity.
* **`satellite`**: Boolean flag to include a moon.
* **`mass_sat`**: Satellite mass [kg].
* **`semimajoraxis_sat`**: Satellite orbital radius [m].

---

### Tidal Model Parameters (`[orbit.obliqua]`)
Controls the tidal response model.
* **`store_3D`**: Boolean flag to store 3D tidal response data.
* **`enforce_ec`**: Boolean flag to enforce energy conservation in tidal response calculations.
* **`optimize_scales`**: Boolean flag to optimize scaling factors for numerical stability.
* **`solid_shell`**: Boolean flag to add an infinitesimal solid shell around the core to couple y2 and y4 in fluid mantles.

#### Rheology and Viscosity
* **`min_frac`**, **`max_frac`**: Minimum segment fraction of total mantle before it is considered.
* **`visc_l`**, **`visc_lus`**: Liquid and liquidus viscosities.
* **`visc_s`**, **`visc_sus`**: Solid and solidus viscosities.

#### Spectral and Forcing Parameters
* **`n`**: Radial dependence exponent in $(r/a)^n$.
* **`m`**: Tidal harmonic (e.g., $m=2$ for semidiurnal tides).
* **`spectrum`**: Frequency sampling strategy (`"full"` or `"adaptive"`).
* **`N_sigma`**: Number of sampled forcing frequencies.
* **`p_min`**, **`p_max`**: Period range ($\log_{10}$ kyr).
* **`s_min`**, **`s_max`**: Fourier mode range.

#### Material Model
* **`material_mu`**: Rheological model for shear modulus (`"andrade"`, `"maxwell"`, or `"elastic"`).
* **`material_k`**: Rheological model for bulk modulus (`"andrade"`, `"maxwell"`, or `"elastic"`).
* **`alpha`**: Andrade power-law exponent.

---

### Interior Structure Models

#### Solid Model (`module_solid`)
* **`solid0d`**: Homogeneous solid approximation.
* **`solid1d`**: Radially resolved structure.
* **`solid1d-relax`**: Relaxation-based solver (more stable).
* **`solid1d-mush`**: Includes partially molten regions.
* **`solid1d-mush-relax`**: Relaxation-based solver (more stable), includes partially molten regions.

#### Mushy Layer (`module_mushy`)
* **`none`**: No explicit treatment.
* **`interp`**: Smooth transition between solid and liquid.

#### Fluid Model (`module_fluid`)
* **`none`**: No fluid layer.
* **`fluid0d`**: Bulk fluid approximation.
* **`fluid1d`**: Radially structured fluid with heating.

---

### Solid Interior Parameters (`[orbit.obliqua.solid]`)
* **`ncalc`**: Number of radial layers (shooting method).
* **`dr_min`**, **`dr_max`**: Grid spacing for relaxation solver [m].
* **`core`**: Core boundary condition (`"liquid"`, `"solid"`, `"inertial"`).
* **`bulk_l`**: Liquid bulk modulus [Pa].
* **`dbulk_power`**: Drained bulk modulus powerlaw scaling exponent.
* **`porosity_thresh`**: Threshold below which no mush is formed.

---

### Fluid Parameters (`[orbit.obliqua.fluid]`)
* **`sigma_R`**: Rayleigh drag at the interface.
* **`sigma_R_inf`**: Drag in the bulk fluid.
* **`sigma_R_prf`**: Vertical drag profile.
* **`H_R`**: Scale height [m].
* **`efficiency`**: Drag efficiency factor.

---

### Mushy Layer Parameters (`[orbit.obliqua.mushy]`)
* **`b_width`**, **`t_width`**: Width of dissipation peaks as fraction of layer thickness.

---

### Planetary Structure (`[struct]`)
Defines bulk planetary properties:
* **`core_density`**: Core density [kg m$^{-3}$].
* **`core_shear`**: Core shear [Pa s].
* **`core_bulk`**: Core bulk modulus [Pa ].

---

### Interior Energetics (`[interior_energetics]`)
Defines interior energy transfer dynamics related properties:
* **`grain_size`**: Grain size [m].