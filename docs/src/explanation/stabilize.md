
### Stabilization techniques for numerical stability
In some cases, the numerical solution of the tidal response can become unstable due to the wide range of physical parameters involved, especially when dealing with fluid-like layers in the solid formalism, that is low-viscosity layers and low forcing frequencies. To mitigate this, we have implemented several stabilization techniques that can be enabled in the configuration file. These techniques are designed to improve the convergence and stability of the numerical solution without significantly affecting the physical accuracy of the results.

### 'enforce_ec' (Enforce Energy Conservation)
When enabled, this option ensures that the computed tidal response adheres to the principle of energy conservation. This is particularly important in the low forcing frequency limit where the tidal response can become highly sensitive to numerical errors. 

The derived heating profile is normalized to ensure that the total energy dissipated matches the expected energy input from the k Lovenumber in a particular segment. This is done by computing the ratio of the expected energy dissipation to the computed energy dissipation and applying this ratio to the heating profile. This order is chosen since the Lovenumbers are generally more stable than the heating profile, which can diverge as y2 and y4 decouple from the rest of the system.

### 'optimize_scales' (Optimize Scaling Factors)
This option allows the model to automatically adjust scaling factors used in the numerical solution to improve stability. The scaling factors are applied to the governing equations to reduce the condition number of the system, which can help prevent numerical instabilities. When this option is enabled, the model will analyze the physical parameters and adjust the scaling factors accordingly.

### 'solid_shell' (Add Infinitesimal Solid Shell)
This option adds an infinitesimal solid shell around the core to couple y2 and y4 in fluid mantles. This is particularly useful in cases where the mantle is mush or fluid-like, as it helps to stabilize the numerical solution by providing a condition that strongly couples the y2 and y4 functions. This generally prevents the decoupling of these functions, which can lead to numerical instabilities in mainly the shear component of the heating profile. The solid shell is infinitesimal and has no complex shear nor bulk modulus, meaning it does not significantly affect the physical properties of the core or mantle, nor does it introduce additional heating in the shell.

