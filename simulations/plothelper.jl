using Plots

mf = font("Roboto-Regular", 10)
lf = font("Roboto-Regular", 8)
default(
    xtickfont = mf,
    ytickfont = mf,
    guidefont = mf,
    legendfont = lf,
    framestyle = :box,
    minorgrid = :false,
    grid = false,
    color = :black
)
