# TunableAerosolCharger

CAD files for the charge strip holder, code for the numerical model, and raw data for the manuscript

> Petters, M.D., Han, S., and Mahant, S.: Design and Characterization of a Modular Tunable Ring Electrode Aerosol Charge Conditioner, J. Aerosol Sci., submitted.

All source code is made available via the MIT License.

## CAD Files

CAD files for printing the holder for the charge strips are in the `cad/` directory. CAD files are released under the CERN Open Hardware Licence Version 2 - Permissive. A copy of the license is provided in this repository.

## Source Code For Numerical Simulations

Code for numerical simulation of the ion field, as well as helper code, are in the `src/` folder. This code is usually accessed through a package interface and can be used for understanding of the numerical implementation or for changing the governing equations used in the simulations.

## Running The Model

Example code for running simulations and extracting solutions is provided in the `simulation/` folder. That code will read an input file specifying parameters and run the code in the `src/` folder. Please refer to the local `README.md` for more information. 

## Data Files

Data files underlying the experiment are in the `data` folder. Each data file has relevant metadata encode into the file as "composition_Size_voltage_polarity.csv". The data correspond to the raw output written by the experimental system.
