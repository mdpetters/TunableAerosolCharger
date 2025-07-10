using TunableAerosolCharger
using YAML
using JLD2

function (@main)(ARGS)
    inputfile = ARGS[1]
    input = YAML.load_file(inputfile)
    out = input["outputfile"]["name"]
    R = input["geometry"]["R"]
    L = input["geometry"]["L"]

    zs, rs, cnewp, cnewm, Ez, Er, Psi = TunableAerosolCharger.simulation(input)

    ii = cnewp .< 0.0
    cnewp[ii] .= 0.0
    jj = cnewm .< 0.0
    cnewm[jj] .= 0.0
    z = zs[3:end-4]
    r = rs[2:end-1]
    Ezs = Ez[3:end-2, 3:end-4]
    Ers = Er[3:end-2, 3:end-4]
    cp = cnewp[3:end-2, 3:end-4] .* 1e6
    cm = cnewm[3:end-2, 3:end-4] .* 1e6
    Psi = Psi[3:end-2, 3:end-4]

    jldsave(out; z = z, r = r, Ezs = Ezs, Ers = Ers, cp = cp, cm = cm, R = R, L = L, Ψ = Psi)
end
