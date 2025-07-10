module TunableAerosolCharger

export attachmentCoefficients

using Random
using QuadGK
using Interpolations
using YAML
using OrdinaryDiffEq
using SteadyStateDiffEq
using Roots
using IntervalRootFinding
using IntervalArithmetic.Symbols
using CSV
using DataFrames
using ForwardDiff
using Memoize
using SpecialFunctions

struct einzelLens
    V1::Float64
    V2::Float64
    D::Float64
    S::Float64
end

str = Base.find_package("TunableAerosolCharger")
p = split(str, "TunableAerosolCharger.jl")[1]
const df⁻ = CSV.read(p*"S1.txt", DataFrame)
const df⁺ = CSV.read(p*"S2.txt", DataFrame)
const df⁰ = CSV.read(p*"S3.txt", DataFrame)

include("fluidFlow.jl")
include("chargingTheory.jl")
include("electricField.jl")
include("dmaRoutines.jl")

end
