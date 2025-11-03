# Getting started
This page outlines requirements and installation steps for the code. Currently,
GNU/Linux and MacOS (including ARM) are supported.

## Software requirements

!!! warning
    Do not install Julia using your system package manager. Install only from julialang.org as below.

## Installation
Follow the steps below in order to setup the code.
1. Install Julia's package manager: `curl -fsSL https://install.julialang.org | sh`
2. Switch to Julia 1.11: `juliaup add 1.11 && juliaup default 1.11`
3. Download Love.jl: `git clone https://github.com/FormingWorlds/Love.jl.git`
4. Change directory: `cd Love.jl`
5. Open Julia package manager: `julia`, then press `]`
6. Install: `add https://github.com/FormingWorlds/Love.jl.git#julia`

## Testing
If you want to run the tests manually, simply use the script in the `test/` folder...
```bash
julia test/runtests.jl
```
This will print information on whether tests passed or failed.

## Using the code
See [Using the model](@ref) for information on using the code.
See [Troubleshooting](@ref) for troubleshooting advice.


