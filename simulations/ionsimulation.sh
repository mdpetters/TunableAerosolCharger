#!/bin/bash

julia simulation.jl input/ionseries1.yml 2>&1 | tee log/ionseries1.txt &
julia simulation.jl input/ionseries2.yml 2>&1 | tee log/ionseries2.txt &
julia simulation.jl input/ionseries3.yml 2>&1 | tee log/ionseries3.txt &
julia simulation.jl input/ionseries4.yml 2>&1 | tee log/ionseries4.txt &
julia simulation.jl input/ionseries5.yml 2>&1 | tee log/ionseries5.txt &
julia simulation.jl input/ionseries6.yml 2>&1 | tee log/ionseries6.txt &
julia simulation.jl input/ionseries7.yml 2>&1 | tee log/ionseries7.txt &
julia simulation.jl input/ionseries8.yml 2>&1 | tee log/ionseries8.txt &
julia simulation.jl input/productionseries1.yml 2>&1 | tee log/productionseries1.txt &
julia simulation.jl input/productionseries2.yml 2>&1 | tee log/productionseries2.txt &
julia simulation.jl input/productionseries3.yml 2>&1 | tee log/productionseries3.txt &
julia simulation.jl input/productionseries4.yml 2>&1 | tee log/productionseries4.txt &
julia simulation.jl input/productionseries5.yml 2>&1 | tee log/productionseries5.txt &
julia simulation.jl input/productionseries6.yml 2>&1 | tee log/productionseries6.txt &
julia simulation.jl input/productionseries7.yml 2>&1 | tee log/productionseries7.txt &
julia simulation.jl input/productionseries8.yml 2>&1 | tee log/productionseries8.txt &
julia simulation.jl input/productionseries9.yml 2>&1 | tee log/productionseries9.txt &
julia simulation.jl input/productionseries10.yml 2>&1 | tee log/productionseries10.txt &

