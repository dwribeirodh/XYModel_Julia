# XYModel_Julia
This package applies the Metropolis algorithm to generate thermodynamic data and equilibrium configurations of the [1D XY model](https://en.wikipedia.org/wiki/Classical_XY_model) and the [2D XY model](https://en.wikipedia.org/wiki/Classical_XY_model#Two_dimensions). Magnetization, internal energy, and specific heat can be computed by running a Monte Carlo simulation with adjustable parameters. Entropy data is computed using the [sweetsourcod](https://github.com/martiniani-lab/sweetsourcod) repo. Make sure to have sweetsourcod functional before trying to compute entropies.

NOTE: this repo is not yet a finished project, it is part of an ongoing research project.
I have some work to do on making the simulation more accessible to other users.

Some functionalities I'm currently testing:
- *Typewriter spin selection (deterministic)*
- *gif generation of simulations*
- *optimizing the lossless compression code*
- *optimizing XY1D_nonper.jl script (currently a little slow)*
- *merging periodic and non-periodic BCs scripts*


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
I have included here some plots generated with the scripts ```XY1D_nonper.jl``` (runs the Monte Carlo simulation). Important simulation parameters:

* epoch=10000000.0
* freq=1000
* L=1000

- *Internal Energy*

![energy_img](https://github.com/dwribeirodh/XYModel_Julia/blob/main/Simulation_Results/2021-08-04/plots/xy_energy_1.0e7_1000.png)

- *Specific Heat Capacity*

![cv_img](https://github.com/dwribeirodh/XYModel_Julia/blob/main/Simulation_Results/2021-08-04/plots/xy_cv_1.0e7_1000.png)


## License

This section will contain info on license
