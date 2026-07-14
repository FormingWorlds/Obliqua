```@meta
CollapsedDocStrings = true
```

### Tutorials (4)

# Plotting 

Once the code has finished running, it will generate and output file that the user can visualize using the python based `visualize_output.ipynb` notebook. This notebook is located in the `examples` folder.

The Obliqua output file will be placed in the `out` folder, unless otherwise specified in the configuration file (default is `out`). The output file will be named according to the timestep of the simulation specified in the configuration file, (default is `0.0_obliqua.nc`).

When you have identified the relevant output file path, you can use the `visualize_output.ipynb` notebook to visualize the output. Simply open the notebook and update the following lines to point to the correct output file path.

```python
out_dir      = "/path/to/obliqua/out/"
nc_file_path = os.path.join(out_dir, "0.0_obliqua.nc")
```

The notebook will then load the output file and automatically produce several plots of the results. These plots will be stored in the `out_dir` directory. 

The user can also use the `Obliqua.plotting` module to generate plots on-the-fly. This is for example useful when troubleshooting. This includes the following functions:
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