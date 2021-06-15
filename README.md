# XYModel_Julia
This package applies the Metropolis algorithm to generate thermodynamic data and equilibrium configurations of the [1D XY model](https://en.wikipedia.org/wiki/Classical_XY_model). Magnetization, internal energy, and specific heat can be computed by running a Monte Carlo simulation with adjustable parameters. Entropy data is computed using the [sweetsourcod](https://github.com/martiniani-lab/sweetsourcod) repo. Make sure to have sweetsourcod functional before trying to compute entropies.

## Installation

To install this package, clone the repository as follows:

```
git clone https://github.com/dwribeirodh/XYModel_Julia.git
```

Some Julia packages are needed to succesfully use the scrips contained in the package. To install the packages, run the following script in the julia REPL.

```
julia> include("/path/to/XYModel_Julia/setup.jl")
```

## Simulation Parameters

To choose different simulation parameters, open the ```config.txt``` file and read the instructions. You can modify the parameters after the "=" sign under each simulation parameter.

## License

This section will contain info on license
