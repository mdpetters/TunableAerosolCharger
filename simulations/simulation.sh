#!/bin/bash

julia simulation.jl input/vseries1.yml 2>&1 | tee log/vseries1.txt &
julia simulation.jl input/vseries2.yml 2>&1 | tee log/vseries2.txt &
julia simulation.jl input/vseries3.yml 2>&1 | tee log/vseries3.txt &
julia simulation.jl input/vseries4.yml 2>&1 | tee log/vseries4.txt &
julia simulation.jl input/vseries5.yml 2>&1 | tee log/vseries5.txt &
julia simulation.jl input/vseries6.yml 2>&1 | tee log/vseries6.txt &
julia simulation.jl input/vseries7.yml 2>&1 | tee log/vseries7.txt &
julia simulation.jl input/vseries8.yml 2>&1 | tee log/vseries8.txt &
julia simulation.jl input/vseries9.yml 2>&1 | tee log/vseries9.txt &
julia simulation.jl input/vseries10.yml 2>&1 | tee log/vseries10.txt &
julia simulation.jl input/vseries11.yml 2>&1 | tee log/vseries11.txt &
julia simulation.jl input/vseries12.yml 2>&1 | tee log/vseries12.txt &
julia simulation.jl input/vseries13.yml 2>&1 | tee log/vseries13.txt &
julia simulation.jl input/vseries14.yml 2>&1 | tee log/vseries14.txt &
julia simulation.jl input/vseries15.yml 2>&1 | tee log/vseries15.txt &
julia simulation.jl input/vseries16.yml 2>&1 | tee log/vseries16.txt &
julia simulation.jl input/vseries17.yml 2>&1 | tee log/vseries17.txt &
julia simulation.jl input/vseries18.yml 2>&1 | tee log/vseries18.txt &
julia simulation.jl input/vseries19.yml 2>&1 | tee log/vseries19.txt &
