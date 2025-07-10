using Interpolations
using Random
using TunableAerosolCharger

d = load("output/vseries9.jld2")
config = YAML.load_file("input/vseries9.yml")

z = d["z"]
r = d["r"]
Ezs = d["Ezs"]
Ers = d["Ers"]
cp = d["cp"]
cm = d["cm"]
R = d["R"]
L = d["L"]

# quick evaluation of Ez and Er using Interpolation
fEz = extrapolate(interpolate((r, z), Ezs, Gridded(Linear())), Flat())
fEr = extrapolate(interpolate((r, z), Ers, Gridded(Linear())), Flat())
fcp = extrapolate(interpolate((r, z), cp, Gridded(Linear())), Flat())
fcm = extrapolate(interpolate((r, z), cm, Gridded(Linear())), Flat())

function ur(r, D, Q)
    R = D / 2.0                # diameter to radius
    q = 1000.0 * Q / 60.0      # L/min to cm3/s
    ub = 4 * q / (pi * D^2.0)  # mean velocity
    u = 2.0 * ub * (1.0 - (r / R)^2.0)
    return u
end

function dzdt(t, r, z, n, B)
    e = 1.60217663e-19
    R = 1.23 * 2.54 / 2.0
    Zp = n * e * B
    dzdt = ur(r, 2 * R, 0.45) - fEz(r, z) * Zp # convert to cm s-1
    return dzdt
end

function drdt(t, r, z, n, B)
    e = 1.60217663e-19
    R = 1.23 * 2.54 / 2.0
    Zp = n * e * B
    drdt = -fEr(r, z) * Zp  # convert to cm s-1
    return drdt
end

function trajectory(r0, z0, nc0, t, h, Dp)
    t, p = 295.15, 1e5
    qsa, qsh = 1.66e-5, 8.3e-5
    r₁, r₂, l = 9.37e-3, 1.961e-2, 0.44369
    Λ = TunableAerosolCharger.DMAconfig(t, p, qsa, qsh, r₁, r₂, l, 0.0, :-, 6, :cylindrical)

    myZ = TunableAerosolCharger.dtoz(Λ, Dp) * 100 * 100 # [cm2 V-1 s-1]
    myB = myZ / 1.60217663e-19
    r = r0
    z = z0
    t0 = 0.0
    m = Int(round((t - t0) / h))
    n = 0
    nc = nc0
    sol = Vector[]
    for i = 1:m
        # z position
        k1 = h * dzdt(t0, r, z, n, myB)
        k2 = h * dzdt(t0 + 0.5 * h, r, z + 0.5 * k1, n, myB)
        k3 = h * dzdt(t0 + 0.5 * h, r, z + 0.5 * k2, n, myB)
        k4 = h * dzdt(t0 + h, r, z + k3, n, myB)

        # r position
        ka = h * drdt(t0, r, z, n, myB)
        kb = h * drdt(t0 + 0.5 * h, r + 0.5 * ka, z, n, myB)
        kc = h * drdt(t0 + 0.5 * h, r + 0.5 * kb, z, n, myB)
        kd = h * drdt(t0 + h, r + kc, z, n, myB)
        z = z + (1.0 / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
        r = r + (1.0 / 6.0) * (ka + 2.0 * kb + 2.0 * kc + kd)

        if (z > 7.6)
            n = 1
        end
        if (z > 20) || (r > R) || (r < -R)
            break
        end
        t0 = t0 + h
        if i % 100 == 0
            push!(sol, [t0, z, r, n])
        end
    end
    return hcat(sol...)
end

config = YAML.load_file("ionProperties.yml")
βd⁺, βd⁻, βe⁺, βe⁻ = TunableAerosolCharger.attachmentCoefficients(100e-9, config)

h = 0.0001
n0 = 0
z0 = 0.6
Dp = 20e-9
sol1 = trajectory(1.25, z0, 0, 40.0, h, Dp)
sol2 = trajectory(1, z0, 0, 40.0, h, Dp)
sol3 = trajectory(0.75, z0, 0, 40.0, h, Dp)
sol4 = trajectory(0.5, z0, 0, 40.0, h, Dp)
sol5 = trajectory(0.25, z0, 0, 40.0, h, Dp)
sol6 = trajectory(0.0, z0, 0, 40.0, h, Dp)
sol7 = trajectory(-0.25, z0, 0, 40.0, h, Dp)
sol8 = trajectory(-0.5, z0, 0, 40.0, h, Dp)
sol9 = trajectory(-0.75, z0, 0, 40.0, h, Dp)
sol10 = trajectory(-1.0, z0, 0, 40.0, h, Dp)
sol11 = trajectory(-1.25, z0, 0, 40.0, h, Dp)

p = heatmap(
    z,
    r,
    (cp .- cm) ./ 1e12,
    cmap = :diverging_bwr_20_95_c54_n256,
    ylabel = "r (cm)",
    xlabel = "z (cm)",
    xlim = [0, L],
    clim = (-3, 3),
    minorticks = true,
)

plot!(sol1[2, :], sol1[3, :], label = :none)
plot!(sol2[2, :], sol2[3, :], label = :none)
plot!(sol3[2, :], sol3[3, :], label = :none)
plot!(sol4[2, :], sol4[3, :], label = :none)
plot!(sol5[2, :], sol5[3, :], label = :none)
plot!(sol6[2, :], sol6[3, :], label = :none)
plot!(sol7[2, :], sol7[3, :], label = :none)
plot!(sol8[2, :], sol8[3, :], label = :none)
plot!(sol9[2, :], sol9[3, :], label = :none)
plot!(sol10[2, :], sol10[3, :], label = :none)
plot!(sol11[2, :], sol11[3, :], label = :none)

Dp = 50e-9
sol1 = trajectory(1.25, z0, 0, 40.0, h, Dp)
sol2 = trajectory(1, z0, 0, 40.0, h, Dp)
sol3 = trajectory(0.75, z0, 0, 40.0, h, Dp)
sol4 = trajectory(0.5, z0, 0, 40.0, h, Dp)
sol5 = trajectory(0.25, z0, 0, 40.0, h, Dp)
sol6 = trajectory(0.0, z0, 0, 40.0, h, Dp)
sol7 = trajectory(-0.25, z0, 0, 40.0, h, Dp)
sol8 = trajectory(-0.5, z0, 0, 40.0, h, Dp)
sol9 = trajectory(-0.75, z0, 0, 40.0, h, Dp)
sol10 = trajectory(-1.0, z0, 0, 40.0, h, Dp)
sol11 = trajectory(-1.25, z0, 0, 40.0, h, Dp)

plot!(sol1[2, :], sol1[3, :], color = :darkgray, ls = :dash, label = :none)
plot!(sol2[2, :], sol2[3, :], color = :darkgray, ls = :dash, label = :none)
plot!(sol3[2, :], sol3[3, :], color = :darkgray, ls = :dash, label = :none)
plot!(sol4[2, :], sol4[3, :], color = :darkgray, ls = :dash, label = :none)
plot!(sol5[2, :], sol5[3, :], color = :darkgray, ls = :dash, label = :none)
plot!(sol6[2, :], sol6[3, :], color = :darkgray, ls = :dash, label = :none)
plot!(sol7[2, :], sol7[3, :], color = :darkgray, ls = :dash, label = :none)
plot!(sol8[2, :], sol8[3, :], color = :darkgray, ls = :dash, label = :none)
plot!(sol9[2, :], sol9[3, :], color = :darkgray, ls = :dash, label = :none)
plot!(sol10[2, :], sol10[3, :], color = :darkgray, ls = :dash, label = :none)
plot!(sol11[2, :], sol11[3, :], color = :darkgray, ls = :dash, label = :none)

p = plot(
    p,
    layout = grid(1, 1),
    size = (500, 250),
    dpi = 300,
    top_margin = 0 * Plots.PlotMeasures.px,
    left_margin = 20 * Plots.PlotMeasures.px,
    right_margin = 10 * Plots.PlotMeasures.px,
    bottom_margin = 10 * Plots.PlotMeasures.px,
)

savefig(p, "f03.png")
FigureServer.pushhtml1("f03.png")
