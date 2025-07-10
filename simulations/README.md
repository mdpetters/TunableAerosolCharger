# Ion-Field Simulations

## Julia Language

All code is written in the Julia language. You can download the current version of the Julia language here: [https://julialang.org/install/](https://julialang.org/install/). The code was developed and tested using Julia version 1.11.5.

The following packages need to be installed before the code can be used.

```julia
julia> using Pkg
julia> Pkg.add(["YAML", "JLD2", "Interpolations", "OrdinaryDiffEq", "Plots", "Roots", "QuadGK", "DataFrames", "CSV", "Random"])
julia> Pkg.add(["https://github.com/mdpetters/TunableAerosolCharger.git"])
```

## Source Code

The source code for the simulations is in the `src/` directory. These files contain the numerical implementation of to solve the governing equations in the manuscript.

## Simulations 

This directory contains code to run simulations. Simulations are initialized using a YAML file which determines the input parameters. Use `basecase.yml` as an example input file. The inputs correspond to the information presented in Table 1 in the manuscript.

### Single Simulation

To run a single simulation from a shell script call

```bash
$ julia simulation.jl input/vseries1.yml
```

This starts the simulation with the input file `vseries1.yml` located in the `input` directory. The output file and directory is defined in the YAML input file. Output is written in `JLD2` file format, which is a hierarchical file format in the Julia language. Simulations may take several hours to complete. After about 1 min, the code will print out the first time step to the terminal. To run simulations in the background, this command can be used. All output is printed to a log file in the log directory.

```bash
$ julia simulation.jl input/vseries1.yml 2>&1 | tee log/vseries1.txt &
```

### Batch simulations

If multiple CPUs are available, it is possible to run simulations in parallel. The script `setup.jl` generates input files varying a single parameter (e.g. voltage or flow). Generate a bash script and run the script to run the model. Edit the file to not exceed the number of physical cores on the machine.

```bash
$ bash simulation.sh  >/dev/null 2>&1
```

# Particle Charging Calculations

The Julia script `chargedistributions.jl` computes the charge distribution at the exit of the charger. The output is written to `csv` files. The `csv` file contains metadata (voltage, diameter, final coordinate, Ls1, and flow rate) as well as charge fractions for charges -30:1:30. The script `chargedistributions.jl` processes charges for a large set of simulations and diameters. Edit the script to change inputs and metadata as needed. 

The `csv` files can also be loaded to further evaluate the charge conditioner output.

# Reading the JLD files

Example scripts generating Figure 2 (`f02.jl`) and Figure 3 (`f03.jl`)in the manuscript are included. These scripts illustrate how to read the JLD files and extract the ion concentration and electric field strength.
