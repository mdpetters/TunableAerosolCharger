using JLD2
using YAML

include("plothelper.jl")

file = "vseries9"
d = load("output/" * file * ".jld2")
config = YAML.load_file("input/" * file * ".yml")
println(config["ions"])
z = d["z"]
r = d["r"]
Ezs = d["Ezs"]
Ers = d["Ers"]
Psi = d["Ψ"]
cp = d["cp"] ./ 1e12
cm = d["cm"] ./ 1e12
R = d["R"]
L = d["L"]

cmi = minimum(cp)
cmi = -3.5
pa = heatmap(
    z,
    r,
    -cm,
    cmap = :diverging_bwr_20_95_c54_n256,
    xlim = [0, L],
    clim = (cmi, abs(cmi)),
    minorticks = true,
    ylabel = "r (cm)",
)

pb = heatmap(
    z,
    r,
    cp,
    cmap = :diverging_bwr_20_95_c54_n256,
    xlim = [0, L],
    clim = (cmi, abs(cmi)),
    minorticks = true,
    ylabel = "r (cm)",
)

pc = heatmap(
    z,
    r,
    cp .- cm,
    cmap = :diverging_bwr_20_95_c54_n256,
    xlabel = "z (cm)",
    xlim = [0, L],
    clim = (cmi, abs(cmi)),
    minorticks = true,
    ylabel = "r (cm)",
)

p1 = contour(
    z,
    r,
    Psi,
    ylabel = "r (cm)",
    xlabel = "",
    cmap = :diverging_bwr_20_95_c54_n256,
    clims = (0, 100),
    xlim = [0, L],
    minorticks = true,
)

p2 = heatmap(
    z,
    r,
    sqrt.(Ezs .^ 2 .+ Ers .^ 2) / 1000,
    cmap = :diverging_rainbow_bgymr_45_85_c67_n256,
    ylabel = "r (cm)",
    xlim = [0, L],
    minorticks = true,
)

n, m = size(cp)
i = Int(n / 2)
p3 = plot(
    z,
    cp[i, :],
    color = :darkred,
    label = "N₊",
    legend = :outertopright,
    xlabel = "z (cm)",
    ylabel = "N (10¹² m⁻³)",
)
p3 = plot!(z, cm[i, :], label = "N₋             ")
xaxis!(minorticks = 4)

p = plot(
    p1,
    pa,
    p2,
    pb,
    p3,
    pc,
    layout = grid(3, 2),
    size = (1100, 500),
    dpi = 300,
    top_margin = 0 * Plots.PlotMeasures.px,
    left_margin = 20 * Plots.PlotMeasures.px,
    right_margin = 0 * Plots.PlotMeasures.px,
    bottom_margin = 20 * Plots.PlotMeasures.px,
)

savefig(p, "f02.png")
# FigureServer.pushhtml1("f02.png")
