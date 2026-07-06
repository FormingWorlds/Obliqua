
# Getting started

This page outlines requirements and installation steps for the code. Currently,
GNU/Linux and MacOS (including ARM) are supported.

## Software requirements

!!! warning
    Do not install Julia using your system package manager. Install only from julialang.org as below.

## Installation

1. Install Julia's package manager (Mac and Linux): 
   ```sh
   curl -fsSL https://install.julialang.org | sh
   ```
2. Clone the repo
   ```sh
   git clone https://github.com/FormingWorlds/Obliqua.git
   ```
3. Move to the repository directory
   ```sh
   cd Obliqua
   ```
4. Open Julia package manager
   ```sh
   julia
   ```
5. Install `Obliqua`
   ```sh
   julia> ]
   pkg> add .
   ```

## Testing
If you want to run the tests manually, simply use the script in the `test/` folder.
```bash
julia test/runtests.jl
```
This will print information on whether tests passed or failed.
