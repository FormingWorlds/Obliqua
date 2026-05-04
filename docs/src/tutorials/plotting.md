```@meta
CollapsedDocStrings = true
```

### Tutorials (4)

# Plotting 

While this is still a work in progress, `Obliqua` currently includes basic plotting functionality. These include the following functions:
- `Obliqua.plotting.plot_imagk2_spectrum`: Plots the imaginary part of the Love number $k_n$ as a function of the forcing frequency.
- `Obliqua.plotting.plot_imagk2_spectra`: Plots the imaginary part of the Love number $k_n$ as a function of the forcing frequency for each segment individually.
- `Obliqua.plotting.save_heat_profile`: Plots the tidal power profile as a function of the radius.
- `Obliqua.plotting.plot_segment_heating`: Plots the tidal power block as a function of the radius and forcing frequency. 
- `Obliqua.plotting.plot_surface_heating`: Plots the surface tidal power as a function of latitude and longitude.
- `Obliqua.plotting.plot_relaxation_solution`: Plots the $y$-functions as a function of the radius. (Specifically useful for debugging the relaxation models, but can be used for the other models as well).


---
---


```@autodocs
Modules = [Obliqua.plotting]
Order   = [:function, :type]   
```

---
---