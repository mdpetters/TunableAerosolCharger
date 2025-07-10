#!/bin/bash

julia simulation.jl input/qseries1.yml 2>&1 | tee log/qseries1.txt &
julia simulation.jl input/qseries2.yml 2>&1 | tee log/qseries2.txt &
julia simulation.jl input/qseries3.yml 2>&1 | tee log/qseries3.txt &
julia simulation.jl input/qseries4.yml 2>&1 | tee log/qseries4.txt &
julia simulation.jl input/qseries5.yml 2>&1 | tee log/qseries5.txt &
julia simulation.jl input/qseries6.yml 2>&1 | tee log/qseries6.txt &
julia simulation.jl input/qseries7.yml 2>&1 | tee log/qseries7.txt &
julia simulation.jl input/qseries8.yml 2>&1 | tee log/qseries8.txt &
julia simulation.jl input/qseries9.yml 2>&1 | tee log/qseries9.txt &
julia simulation.jl input/qseries10.yml 2>&1 | tee log/qseries9.txt &
julia simulation.jl input/qseries11.yml 2>&1 | tee log/qseries1.txt &
julia simulation.jl input/qseries12.yml 2>&1 | tee log/qseries2.txt &
julia simulation.jl input/qseries13.yml 2>&1 | tee log/qseries3.txt &
julia simulation.jl input/qseries14.yml 2>&1 | tee log/qseries4.txt &
julia simulation.jl input/qseries15.yml 2>&1 | tee log/qseries5.txt &
julia simulation.jl input/qseries16.yml 2>&1 | tee log/qseries6.txt &
julia simulation.jl input/qseries17.yml 2>&1 | tee log/qseries7.txt &
julia simulation.jl input/qseries18.yml 2>&1 | tee log/qseries8.txt &
julia simulation.jl input/qseries19.yml 2>&1 | tee log/qseries9.txt &
julia simulation.jl input/qseries20.yml 2>&1 | tee log/qseries9.txt &
