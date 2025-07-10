# ---------------------------------------------------------------------------------
# fluidFlow.jl
#
# Calculate the flow velocity profile in the tube
#
# Author: Markus Petters (markus.petters@ucr.edu)
# ---------------------------------------------------------------------------------

"""
    function u(r, D, Q)

Inputs:

- r      radial position [cm]
- D      tube diameter [cm]
- Q      volumetric flow rate [L min-1]

Outputs:

- u      flow velocity [cm s-1]
"""
function ur(r, D, Q)
    R = D / 2.0                # diameter to radius
    q = 1000.0 * Q / 60.0      # L/min to cm3/s
    ub = 4 * q / (pi * D^2.0)  # mean velocity
    u = 2.0 * ub * (1.0 - (r / R)^2.0)
    return u
end

"""
    function simulation(input)

Inputs:
- input     Model input parsed from YAML file

Outputs:
- zs        Grid points in flow dimension               [cm]
- rs        Grid points in radial dimension             [cm]
- cp        Positive ion concentration                  [# cm-3]
- cm        Negative ion concentration                  [# cm-3]
- Ez        Electric field strength in flow direction   [V cm-1]
- Er        Electric field strength in radial direction [V cm-1]
- Psi       Potential field                             [V]

References:

- Advection Upwind scheme from Fromm (1968) as summarized in lecture notes by van den Bosch (2021) Numerical Hydrodynamics, http://www.astro.yale.edu/vdbosch/Numerical_Hydrodynamics.pdf
- Diffusion Second order finite differences as described in Numerical Recipes (Press et al. 1992).

Fromm, Jacob E. "A method for reducing dispersion in convective difference schemes." Journal of Computational Physics 3, no. 2 (1968): 176-189.

Press, William H.; Teukolsky, Saul A.; Vetterling, William T.; Flannery, Brian P. (1902). Numerical Recipes: The Art of Scientific Computing.
"""
function simulation(input)
    dt = input["simulation"]["dt"]      # [s]
    tmax = input["simulation"]["tmax"]  # [s]
    n = input["simulation"]["n"]        # [-]
    m = input["simulation"]["m"]        # [-]

    Zpp = input["ions"]["Zpp"]          # [cm2 V-1 s-1]
    Zpm = input["ions"]["Zpm"]          # [cm2 V-1 s-1]
    α = input["ions"]["α"]              # [cm3 s-1]
    ionp = input["ions"]["p"]           # [cm-3 s-1]
    D = input["ions"]["D"]              # [cm2 s-1]

    v1 = input["potential"]["v1"]       # [V]
    v2 = input["potential"]["v2"]       # [V]

    R = input["geometry"]["R"]          # [cm]
    L = input["geometry"]["L"]          # [cm]
    Lb1 = input["geometry"]["Lb1"]      # [cm]
    Lb2 = input["geometry"]["Lb2"]      # [cm]
    Ls1 = input["geometry"]["Ls1"]      # [cm]
    Ls2 = input["geometry"]["Ls2"]      # [cm]
    S = input["geometry"]["S"]          # [cm]

    q = input["flow"]["rate"]           # [L min-1]

    rs = range(-R, stop = R, length = n - 2)
    zs = range(0, stop = L, length = m)
    j1 = findfirst(x -> x > Ls1, zs)
    j2 = findfirst(x -> x > Ls2, zs)
    dz = zs[2] - zs[1]
    c0 = zeros(n, m)
    dr = rs[3] - rs[2]
    r = [rs[1] - dr; rs; rs[end] + dr]
    u = [0.0; ur.(rs, 2 * R, q); 0.0]
    lens1 = einzelLens(v2, v1, 2.0 * R, S)
    lens2 = einzelLens(v1, v2, 2.0 * R, S)
    Psia = zeros(n, m)
    Eza = zeros(n, m)
    Era = zeros(n, m)
    Psi1, Ez1, Er1, _, _ = lensFields(lens1, (zs .- Lb1), rs)
    Psi2, Ez2, Er2, _, _ = lensFields(lens2, (zs .- Lb2), rs)
    Psi = Psi1 .+ Psi2 .- v1
    Ez = Ez1 .+ Ez2
    Er = Er1 .+ Er2
    Psia[2:(end - 1), :] .= Psi
    Eza[2:(end - 1), :] .= Ez
    Era[2:(end - 1), :] .= Er
    Psi = Psia
    Ez = Eza
    Er = Era

    function f!(cnew, c, dt, Zp)
        maxfac = 0.0
        for i = 3:(n - 2)
            for j = 3:(m - 2)
                vzfac = -Zp * Ez[i, j] * dt / dz
                ufac = u[i] * dt / dz + vzfac
                # ufac = vzfac
                (ufac < maxfac) || (maxfac = ufac)
                if ufac > 0
                    cnew[i, j] =
                        cnew[i, j] -
                        0.25 *
                        ufac *
                        (c[i, j + 1] + 3.0 * c[i, j] - 5.0 * c[i, j - 1] + c[i, j - 2]) +
                        0.25 * ufac^2 * (c[i, j + 1] - c[i, j] - c[i, j - 1] + c[i, j - 2])
                else
                    cnew[i, j] =
                        cnew[i, j] +
                        0.25 *
                        ufac *
                        (c[i, j - 1] + 3.0 * c[i, j] - 5.0 * c[i, j + 1] + c[i, j + 2]) +
                        0.25 * ufac^2 * (c[i, j - 1] - c[i, j] - c[i, j + 1] + c[i, j + 2])
                end
                vrfac = -Zp * Er[i, j] * dt / (r[i] - r[i - 1])
                (vrfac < maxfac) || (maxfac = vrfac)
                if (vrfac > 0)
                    cnew[i, j] =
                        cnew[i, j] -
                        0.25 *
                        vrfac *
                        (c[i + 1, j] + 3.0 * c[i, j] - 5.0 * c[i - 1, j] + c[i - 2, j]) +
                        0.25 * vrfac^2 * (c[i + 1, j] - c[i, j] - c[i - 1, j] + c[i - 2, j])
                else
                    cnew[i, j] =
                        cnew[i, j] +
                        0.25 *
                        vrfac *
                        (c[i - 1, j] + 3.0 * c[i, j] - 5.0 * c[i + 1, j] + c[i + 2, j]) +
                        0.25 * vrfac^2 * (c[i - 1, j] - c[i, j] - c[i + 1, j] + c[i + 2, j])
                end
                # Diffusion in z direction
                zfac = D * dt / dz^2
                (zfac < maxfac) || (maxfac = zfac)
                cnew[i, j] = cnew[i, j] + zfac * (c[i, j + 1] - 2 * c[i, j] + c[i, j - 1])

                # Diffusion in r direction
                g(i) = r[i] * (c[i + 1, j] - c[i - 1, j]) / (r[i + 1] - r[i - 1])
                cnew[i, j] = cnew[i, j] + D * dt / r[i] * (g(i + 1) - g(i - 1)) / (2.0 * dr)
                if (j >= j1) & (j <= j2)
                    cnew[i, j] += ionp * dt
                end
            end
        end
        return maxfac
    end

    t = 0.0
    cnewp = deepcopy(c0)
    cnewm = deepcopy(c0)
    i = 1
    facp = 0.0
    facm = 0.0
    count = 0
    while (t < tmax)
        if (isnan.(cnewp) |> sum) > 0
            println((isnan.(cnewp) |> sum))
            break
        end
        if (count % 1000) == 0
            println(t)
        end
        cp = deepcopy(cnewp)
        cm = deepcopy(cnewm)
        facp = f!(cnewp, cp, dt, Zpp)
        facm = f!(cnewm, cm, dt, Zpm)
        recomb = cnewp .* cnewm .* α .* dt
        ii = (cnewp .< 0) .| (cnewm .< 0)
        recomb[ii] .= 0
        for i = 1:n
            for j = 1:m
                if (recomb[i, j] < cnewp[i, j]) & (cnewp[i, j] > 0.0) & (recomb[i, j] > 0)
                    cnewp[i, j] = cnewp[i, j] - recomb[i, j]
                end
                if (recomb[i, j] < cnewm[i, j]) & (cnewm[i, j] > 0.0) & (recomb[i, j] > 0)
                    cnewm[i, j] = cnewm[i, j] - recomb[i, j]
                end
            end
        end
        t = t + dt
        count = count + 1
    end
    println(facp)
    println(t)
    return zs, rs, cnewp, cnewm, Ez, Er, Psi
end


