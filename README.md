# XYModel_Julia
This package applies the Metropolis algorithm to generate thermodynamic data and equilibrium configurations of the [1D XY model](https://en.wikipedia.org/wiki/Classical_XY_model). Magnetization, internal energy, and specific heat can be computed by running a Monte Carlo simulation with adjustable parameters. Entropy data is computed using the [sweetsourcod](https://github.com/martiniani-lab/sweetsourcod) repo. Make sure to have sweetsourcod functional before trying to compute entropies.

## Installation

To install this package, clone the repository as follows:

```
git clone https://github.com/dwribeirodh/XYModel_Julia.git
```

Some Julia packages are needed to successfully use the scrips contained in the package. To install the packages, run the following command in the julia REPL.

```
julia> include("/path/to/XYModel_Julia/setup.jl")
```

## Simulation Parameters

To choose different simulation parameters, open the ```config.txt``` file and read the instructions. You can modify the parameters after the "=" sign under each simulation parameter.

## Simulation Results
I have included here some plots generated with the scripts ```XY1D_nonper.jl``` (runs the Monte Carlo simulation),
and ```XY1D_nonper_entropy.jl``` (lossless compression and entropy computation).

- *Internal Energy*

![energy_img](https://github.com/dwribeirodh/XYModel_Julia/blob/main/Simulation_Results/2021-06-28/plots/xy_energy_5.0e7_100000.png)

- *Specific Heat Capacity*

![cv_img](https://github.com/dwribeirodh/XYModel_Julia/blob/main/Simulation_Results/2021-06-28/plots/xy_cv_5.0e7_100000.png)

- *Entropy*

![entropy_img](https://github.com/dwribeirodh/XYModel_Julia/blob/main/Simulation_Results/2021-06-28/plots/1d_xy_entropy.png)

The simulation parameters used to obtain these results were:
      * L = 100000
      * freq = 100
      * epoch = 50000000.0

## License

This section will contain info on license
