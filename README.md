

<a href="https://opensource.org/license/mit">
   <img src="https://img.shields.io/badge/License-MIT-blue.svg">
</a>
<a href="https://arxiv.org/abs/2507.11266">
  <img src="https://img.shields.io/badge/arXiv-2507.11266-b31b1b">
</a>

## Obliqua (Tidal heating model)

A Julia package to calculate the tidal deformation (i.e., tidal Love numbers) of solid, partially-solid, and liquid planetary mantles

Forked from the [original repository](https://github.com/hamishHay/Love.jl) of Hamish Hay. Distributed under the MIT License.

### Documentation
https://proteus-framework.org/Obliqua

## Contributors

| Name                    | Email address                         |
| -                       | -                                     |
| Marijn van Dijk         | m.r.van.dijk.3[at]student.rug.nl      |
| Tim Lichtenberg         | tim.lichtenberg[at]rug.nl             |
| Harrison Nicholls       | harrison.nicholls[at]physics.ox.ac.uk |
| Mohammad Farhat         | farhat[at]berkeley.edu                |
| Hamish Hay              | hamish[at]tides.rocks                 |


### Repository structure
* `README.md`           - This file
* `LICENSE.txt`         - License for modification, distribution, etc.
* `src/`                - Source files
* `examples/`           - Tools and notebooks showing how to use the code

### Installation

1. Clone the repo
   ```sh
   git clone https://github.com/FormingWorlds/Obliqua.git
   ```
2. Move to the repository directory and start Julia
   ```sh
   cd Obliqua
   julia
   ```
3. Install `Obliqua`
   ```sh
   julia> ]
   pkg> add https://github.com/FormingWorlds/Obliqua.git#julia

### Integrations

This project is being developed in close collaboration with the [PROTEUS](https://github.com/FormingWorlds/PROTEUS) similation framework. Obliqua can be used seamlessly within PROTEUS to model tidal dissipation and Love number evolution in planetary interiors.

For detailed instructions on installing and using Obliqua in combination with PROTEUS, please refer to the [installation guide](https://proteus-framework.org/Obliqua/dev/install/) in the Obliqua documentation.
