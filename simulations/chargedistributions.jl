using TunableAerosolCharger
using YAML
using JLD2
using Interpolations
using OrdinaryDiffEq
using Plots
using Roots
using QuadGK
using DataFrames
using CSV

include("plothelper.jl")

function ur(r, D, Q)
    R = D / 2.0                # diameter to radius
    q = 1000.0 * Q / 60.0      # L/min to cm3/s
    ub = 4 * q / (pi * D^2.0)  # mean velocity
    u = 2.0 * ub * (1.0 - (r / R)^2.0)
    return u
end

function qcyl(r1, r2, R, Q)
    out = quadgk(r -> 2π * r * ur(r, 2 * R, Q), r1, r2, rtol = 1e-9)
    return out[1]
end

function iso(n, R, Q)
    dv = qcyl(0.0, R, R, Q) ./ n

    tr = [0.0]
    for i = 1:(n - 1)
        myr = find_zero(r -> qcyl(tr[i], r, R, Q) - dv, (0.0, R))
        push!(tr, myr)
    end
    tr = [tr; R]
    return tr[1:(end - 1)] .+ 0.5 * (tr[2:end] .- tr[1:(end - 1)])
end

function tunableCharging(r, Dd, z, R, Q, fcp, fcm, zend)
    config = YAML.load_file("ionProperties.yml")

    k = 32
    x = zeros(k)
    Nin = [1.0; x; x]
    zpos(t) = ur(r, 2.0 * R, Q) * t
    n⁺(t) = fcp(r, zpos(t))
    n⁻(t) = fcm(r, zpos(t))

    ts = 0:0.01:30
    out = zend .- map(zpos, ts)
    ~, i = findmin(out .^ 2)
    tend = ts[i]

    βd⁺, βd⁻, βe⁺, βe⁻ = TunableAerosolCharger.attachmentCoefficients(Dd, config)
    param = (βd⁺, βd⁻, n⁺, n⁻, k)
    problem = ODEProblem(TunableAerosolCharger.birthAndDeath, Nin, (0.0, tend), param)
    solution = solve(problem, RK4(), reltol = 1e-4, abstol = 1e-4)
    data = hcat(solution.u...)
    f₋ = data[(2 + k):(end - 1), end]
    f₊ = data[2:k, end]
    f₀ = data[1, end]
    return (f₀, f₊, f₋)
end

function chargeModel(Dd, file, zend)
    d = load("output/" * file * ".jld2")
    config = YAML.load_file("input/" * file * ".yml")

    z = d["z"]
    r = d["r"]
    Ezs = d["Ezs"]
    Ers = d["Ers"]
    cp = d["cp"]
    cm = d["cm"]
    R = d["R"]
    L = d["L"]
    Q = config["flow"]["rate"]
    v = abs(config["potential"]["v1"])

    # quick evaluation of Ez and Er using Interpolation
    fEz = extrapolate(interpolate((r, z), Ezs, Gridded(Linear())), Flat())
    fEr = extrapolate(interpolate((r, z), Ers, Gridded(Linear())), Flat())
    fcp = extrapolate(interpolate((r, z), cp, Gridded(Linear())), Flat())
    fcm = extrapolate(interpolate((r, z), cm, Gridded(Linear())), Flat())

    nt = 100
    rs = iso(nt, R, Q)
    if v > 500
        x = 90
    else
        x = 25
    end
    data = map(1:(nt - x)) do i
        r = rs[i]
        out = tunableCharging(r, Dd, z, R, Q, fcp, fcm, zend)
        f₀ = out[1]
        f₊ = out[2]
        f₋ = out[3]
        (f₀, f₊, f₋)
        [reverse(f₋[1:30]); f₀; f₊[1:30]]
    end

    a = hcat(data...)
    theory = sum(a, dims = 2) ./ (nt - x)
    ns = -30:1:30
    a = DataFrame(
        v = config["potential"]["v1"],
        Dd = Dd,
        zend = zend,
        Ls1 = config["geometry"]["Ls1"],
        q = config["flow"]["rate"],
    )
    b = DataFrame(theory', :auto)
    out = hcat(a, b)
    return out
end

files = [
    "vseries19",
    "vseries18",
    "vseries17",
    "vseries16",
    "vseries15",
    "vseries14",
    "vseries13",
    "vseries12",
    "vseries11",
    "vseries10",
    "vseries9",
    "vseries8",
    "vseries7",
    "vseries6",
    "vseries5",
    "vseries4",
    "vseries3",
    "vseries2",
    "vseries1",
]
Dds = [
    10e-9,
    20e-9,
    30e-9,
    40e-9,
    50e-9,
    60e-9,
    70e-9,
    80e-9,
    90e-9,
    100e-9,
    200e-9,
    300e-9,
    400e-9,
    500e-9,
]

data = map(Dds) do Dd
    intermediate = map(files) do file
        println(Dd, " ", file)
        chargeModel(Dd, file, 18.7)
    end
    vcat(intermediate...)
end
df = vcat(data...)
df |> CSV.write("charging_voltage.csv")

files = [
    "ionseries8",
    "ionseries7",
    "ionseries6",
    "ionseries5",
    "ionseries4",
    "ionseries3",
    "ionseries2",
    "ionseries1",
]

data = map(Dds) do Dd
    intermediate = map(files) do file
        println(Dd, " ", file)
        chargeModel(Dd, file, 18.7)
    end
    vcat(intermediate...)
end
df = vcat(data...)
df |> CSV.write("charge_ions.csv")

files = [
    "productionseries10",
    "productionseries9",
    "productionseries8",
    "productionseries7",
    "productionseries6",
    "productionseries5",
    "productionseries4",
    "productionseries3",
    "productionseries2",
    "productionseries1",
]
data = map(Dds) do Dd
    intermediate = map(files) do file
        println(Dd, " ", file)
        chargeModel(Dd, file, 18.7)
    end
    vcat(intermediate...)
end
df = vcat(data...)
df |> CSV.write("charge_production.csv")

files = [
    "qseries20",
    "qseries19",
    "qseries18",
    "qseries17",
    "qseries16",
    "qseries15",
    "qseries14",
    "qseries13",
    "qseries12",
    "qseries11",
    "qseries10",
    "qseries9",
    "qseries8",
    "qseries7",
    "qseries6",
    "qseries5",
    "qseries4",
    "qseries3",
    "qseries2",
    "qseries1",
]
data = map(Dds) do Dd
    intermediate = map(files) do file
        println(Dd, " ", file)
        chargeModel(Dd, file, 18.7)
    end
    vcat(intermediate...)
end
df = vcat(data...)
df |> CSV.write("charge_flow.csv")
