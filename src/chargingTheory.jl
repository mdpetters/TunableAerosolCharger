# ---------------------------------------------------------------------------------
# chargingTheory.jl
#
# Functions to calculate particle charging
#
# Author: Markus Petters (markus.petters@ucr.edu)
# ---------------------------------------------------------------------------------

"""
    function calculateFieldBeta(E0, Dp, n, zi, t, p, ϵ)

    Purpose:
Computes the ion-particle attachement coefficient for aerosol field charging based on Eq. 8 in He et al. (2015).

    Inputs:
- E0    Electric field strength [V m-1]
- Dp    Particle Diameter [m]
- n     Number of particle charges [-]
- Zi    Ion mobility [m2 V-1 s-1]
- T     Temperature [K]
- p     Pressure [Pa]
- ϵ     Dielectric contant of particle [-]

    Outputs:
- βe    Ion–particle attachment coefficient [m3 s-1]

    Reference
He, K., Lu, J., Ma, X., Ju, Y., Xie, L., Pang, L., Wang, X., & Chen, J. (2015). Effect of Maxwell–Wagner Relaxation on Field Charging of Particles. Aerosol Science and Technology, 49(12), 1210–1221. https://doi.org/10.1080/02786826.2015.1112875
"""
function calculateFieldBeta(E0, Dp, n, Zi, T, p, ϵ)
    e = 1.60217662e-19    # Elementary charge [C]
    kb = 1.380649e-23     # Boltzmann constant [J K-1]
    KE = 9.0 * 1e9        # constant of proportionality [N m2 C-2]
    a = Dp / 2.0          # radius
    ki = Zi * p / 101315.0  # nomenclature

    # Saturation charge Eq. (15.25) Hinds and Zhu 2022
    ns = 3.0 * ϵ / (ϵ + 2.0) * E0 * a^2 / (KE * e)
    if n ≥ ns
        return 0.0
    else
        return 3.0 * ϵ / (ϵ + 2.0) * π * a^2 * ki * E0 * (1.0 - n / ns)^2
    end
end


"""
    function calculateFuchsBeta(Dp, n, Zi, t, p, ϵ)

    Purpose:
Computes the ion-particle attachement coefficient for aerosol field charging based aon the limiting sphere theory (Fuchs 1963).

    Inputs:
- Dp    Particle Diameter [m]
- n     Number of particle charges [-]
- Zi    Ion mobility [m2 V-1 s-1]
- mi    Mass per ion [kg]
- T     Temperature [K]
- p     Pressure [Pa]
- ϵ     Dielectric contant of particle [-]

    References:

Unipolar diffusion charging of aerosol particles in the transition regime. Journal of Aerosol Science, 36(2), 247–265. https://doi.org/10.1016/j.jaerosci.2004.09.002

Fuchs, N.A. On the stationary charge distribution on aerosol particles in a bipolar ionic atmosphere. Geofisica Pura e Applicata 56, 185–193 (1963). https://doi.org/10.1007/BF01993343.

D. Y. H. Pui, S. Fruin & P. H. McMurry (1988) Unipolar Diﬀusion Charging of Ultraﬁne Aerosols, Aerosol Science and Technology, 8:2, 173-187, DOI: 10.1080/02786828808959180.
"""
function calculateFuchsBeta(Dp, n, Zi, mi, T, p, ϵ)
    e = 1.60217662e-19          # Elementary charge [C]
    kb = 1.380649e-23           # Boltzmann constant [J K-1]
    KE = 9.0 * 1e9              # constant of proportionality [N m2 C-2]
    mj = 28.97 / 6.023e23         # mass per molecule of dry air [kg]

    a = Dp / 2                  # Particle radius [m]
    ki = Zi * p / 101315.0      # Pressure dependence of ion mobility (Biskos et al. 2005)
    Di = kb * T * ki / e        # Diffusion coefficient Eq. 10 (Pui et al. 1988) [m2 s-1]
    ci = sqrt(8.0 * kb * T / (π * mi))  # Mean thermal speed Eq 11 (Pui et al. 1988) [m s-1]

    # Mean free path Eq. 12 Pui et al. (1988)
    λ = 8.0 / (3.0 * sqrt(π)) * Di / (1 + 0.016) * sqrt(8.0 / π) / ci / sqrt((mi + mj) / mj) # [m]
    λ = 1.329 * ki / e * sqrt(kb * T * mi * mj / (mi + mj))
    # Limiting sphere radius (Fuchs 1963)
    δ1 = ((1 + (λ / a))^5) / 5
    δ2 = ((1 + (λ^2 / a^2)) * (1 + (λ / a))^3) / 3
    δ3 = 2 / 15 * (1 + λ^2 / a^2)^2.5
    δ = (a^3 / λ^2) * (δ1 - δ2 + δ3)

    # ϕ(r) in Eq. 4 in Biskos et al. 2005
    κ = (ϵ - 1) * e^2 / (ϵ + 1) # image force parameter for particles with ϵ
    ϕ(n, r) = KE * ((n * e^2 / r) - (κ * a^3 / (2 * r^2 * (r^2 - a^2))))

    # prob. of ion entering the limiting-sphere to collide and transfer charge to particle
    b2(r) = r^2 * (1.0 + 2.0 / (3.0 * kb * T) * (ϕ(n, δ) - ϕ(n, r)))
    bm = minimum(b2.(a:(a / 10000.0):δ)) |> sqrt
    γ = bm^2 / δ^2
    if γ ≥ 1.0
        γ = 1.0
    end

    # Integral in denominator of Eq. 11 (Fuchs 1963)
    g(x) = exp(ϕ(n, a / x) / (kb * T))
    Ψ = quadgk(g, 0, a / δ, rtol = 1e-9)[1]

    # Ion-particle attachement coefficient, Eq. 10 (Fuchs 1963)
    β =
        (4.0 * π * a * Di) /
        (((4 * Di * a) / (γ .* ci .* δ^2) * exp(ϕ(n, δ) / (kb * T))) + Ψ)
    βd = (β > 0) ? β : 0.0
    return βd
end

"""
    tryβ(f::Function, val)

Convenience function for try catch block, returning 0 upon failure.
"""
function tryβ(f::Function, val)
    try
        f(val)
    catch
        0.0
    end
end

"""
    attachementCoefficients(Dp::Float64, config::Dict)

    Inputs:
- Dp       Particle diameter [m]
- config   Configuration as dictionary obtained from YAML file

    Outputs:
Lookup function for ion-particle attachment coefficients. n is the number of charges on the particle (-100:1:100), E0 is the electric field strength [V m-1]. The functions have the ion properties, thermodynamic state, and dielectric constant applied. For diffusion, a lookup function is provided for fast retrieval. For field charging the calculateFieldBeta is called.

- βd⁺(n)       positive ions-particle attachment diffusion Fuchs limiting sphere.
- βd⁻(n)       negative ions-particle attachment diffusion Fuchs limiting sphere.
- βe⁺(E0, n)   positive ions-particle attachment field charging.
- βe⁻(E0, n)   negative ions-particle attachment field charging.

After applying the function the attachment coefficient β in [m3 s-1] is returned

    Example:

Obtain coefficients for 100 nm particles

```julia
julia> config = YAML.load_file("configuration.yml")
julia> βd⁺, βd⁻, βe⁺, βe⁻ = attachementCoefficients(100e-9, config)
julia> βd⁻(-2)
5.613089412788073e-13
```
"""
function attachmentCoefficients(Dp::Float64, config::Dict)
    T = config["thermo"]["T"]                        # [K]
    p = config["thermo"]["p"]                        # [Pa]
    Zpp = config["ions"]["Zpp"] * 1e-4               # [m2 V-1 s-1]
    Zpm = config["ions"]["Zpm"] * 1e-4               # [m2 V-1 s-1]
    mip = config["ions"]["mip"] / 1000.0 / 6.023e23  # [kg]
    mim = config["ions"]["mim"] / 1000.0 / 6.023e23  # [kg]
    ϵ = config["particle"]["epsilon"]                # [-]

    # positive ion-particle attachment coefficient
    # mβd⁺(n) = calculateFuchsBeta(Dp, n, Zpp, mip, T, p, ϵ)
    mβd⁺(n) = calculateXerxesBeta(Dp, n, df⁺)
    # negative ion-particle attachment coefficient
    # mβd⁻(n) = calculateFuchsBeta(Dp, -n, abs(Zpm), mim, T, p, ϵ)
    mβd⁻(n) = calculateXerxesBeta(Dp, n, df⁻)

    # Create lookup/convenience functions for β to speed up calculation
    charges = -100:1:100
    mβp = map(n -> tryβ(mβd⁺, n), charges)
    βd⁺ = interpolate((charges,), mβp, Gridded(Linear()))
    mβm = map(n -> tryβ(mβd⁻, n), charges)
    βd⁻ = interpolate((charges,), mβm, Gridded(Linear()))
    βe⁺(E0, n) = calculateFieldBeta(E0, Dp, n, Zpp, T, p, ϵ)
    βe⁻(E0, n) = calculateFieldBeta(E0, Dp, -n, abs(Zpm), T, p, ϵ)

    return βd⁺, βd⁻, βe⁺, βe⁻
end

"""
    birthAndDeath(N, p, t)

    Purpose:
Function to calculate the generation of charged particles for diffusion charging in an ion field. Interfaces with OrdinaryDifferentialEq package to solve for the time evolution. Equations are from Hoppel (1985), Eqs. (3)-(5).

    Inputs:
- p        Parameters
    1. p[1] : ion attachment coefficient function for positively charged particles
    2. p[2] : ion attachment coefficient function for negatively charged particles
    3. p[3] : positive ion concentration [# m-3] as function of time
    4. p[4] : negative ion concentration [# m-3] as function of time
    5. p[5] : upper limit of charges considered
- t        time variable for DiffEq solver
- N        Vector of particle concentrations [N₀; N₊; N₋]
    1. N₀   : neutral particles scalar
    2. N₊   : vector 1..k positively charged particles
    3. N₋   : vector 1..k negatively charged particles


    Outputs:

- dN       Vector of derivatives for DiffEq solver [dN₀; dN₊; dN₋]


    References:
Hoppel, W. A. (1985), Ion-aerosol attachment coefficients, ion depletion, and the charge distribution on aerosols, J. Geophys. Res., 90(D4), 5917–5923, doi:10.1029/JD090iD04p05917.

    Example:

```julia
    config = yaml.load_file("ionproperties.yml")
    k = 12sdasd
    x = zeros(k)
    Nin = [1.0; x; x]
    tspan = (0.0, 5.0)
    Dd = 50e-9
    βd⁺, βd⁻, βe⁺, βe⁻ = tunableaerosolcharger.attachmentcoefficients(Dd, config)
    param = (βd⁺, βd⁻, 1e12, 1e12, k)
    problem = odeproblem(TunableAerosolCharger.birthAndDeath, Nin, tspan, param)
    solution = solve(problem, rk4(), reltol = 1e-8, abstol = 1e-8)
    data = hcat(solution.u...)
    f₀ = data[1, end]
    f₊ = data[2:k, end]
    f₋ = data[2+k:end-1, end]
```
"""
function birthAndDeath(N, p, t)
    β₊ = p[1]
    β₋ = p[2]
    n₊ = p[3](t)
    n₋ = p[4](t)
    i = Int(p[5])
    N₀ = N[1]
    N₊ = N[2:(2 + i - 1)]
    N₋ = N[(2 + i):(2 + i + i - 1)]
    dN₊ = zeros(i)
    dN₋ = zeros(i)

    dN₀ = n₊ * β₊(-1) * N₋[1] + n₋ * β₋(1) * N₊[1] - n₊ * β₊(0) * N₀ - n₋ * β₋(0) * N₀
    dN₊[1] = n₊ * β₊(0) * N₀ - n₊ * β₊(1) * N₊[1] + n₋ * β₋(2) * N₊[2] - n₋ * β₋(1) * N₊[1]
    dN₋[1] =
        n₋ * β₋(0) * N₀ - n₋ * β₋(-1) * N₋[1] + n₊ * β₊(-2) * N₋[2] - n₊ * β₊(-1) * N₋[1]

    for k = 2:(i - 1)
        dN₊[k] =
            n₊ * β₊(k - 1) * N₊[k - 1] - n₊ * β₊(k) * N₊[k] + n₋ * β₋(k + 1) * N₊[k + 1] -
            n₋ * β₋(k) * N₊[k]

        dN₋[k] =
            n₋ * β₋(-(k - 1)) * N₋[k - 1] - n₋ * β₋(-k) * N₋[k] +
            n₊ * β₊(-(k + 1)) * N₋[k + 1] - n₊ * β₊(-k) * N₋[k]
    end

    return [dN₀; dN₊; dN₋]
end

function calculateXerxesBeta(Dp, n, df)
    a = Dp / 2.0
    charges = df[!, 1]
    ii = charges .== n
    k = findmax(ii)[2]
    B = map(i -> df[i, 2:(end - 4)] |> Vector, 1:length(charges))
    n = length(B[1]) - 1
    ex = mapfoldl(i -> B[k][i + 1] * (log10(a))^i, +, 0:n)
    β = exp10(ex)

    return β
end

function calculateXerxesFraction(Dp, n, df)
    a = Dp / 2.0
    charges = df[!, 1]
    ii = charges .== n
    k = findmax(ii)[2]
    B = map(i -> df[i, 2:(end - 4)] |> Vector, 1:length(charges))
    n = length(B[1]) - 1
    ex = mapfoldl(i -> B[k][i + 1] * (log10(a))^i, +, 0:n)
    f = exp10(ex)

    return f
end


"""
   equilibriumDiffusionCharging(Dd, config)

    Purpose:
Compute the equilibrium charging based on the ion attachment coefficients, the birth and death model equations, and balanced ion concentrations.

    Inputs:
- Dp       Particle diameter [m]
- config   Configuration as dictionary obtained from YAML file


    Output:
- f        Vector of charge fractions [f₀; f₊; f₋]
    1. f₀   : neutral particles scalar
    2. f₊   : vector 1..k positively charged particles
    3. f₋   : vector 1..k negatively charged particles

    Example:

```julia
julia>
```
"""
function equilibriumDiffusionCharging(Dd, config)
    k = 15
    x = zeros(k)
    Nin = [1.0; x; x]
    n⁺(t) = 1e12
    n⁻(t) = 1e12
    βd⁺, βd⁻, βe⁺, βe⁻ = TunableAerosolCharger.attachmentCoefficients(Dd, config)
    param = (βd⁺, βd⁻, n⁺, n⁻, k)
    problem = ODEProblem(TunableAerosolCharger.birthAndDeath, Nin, (0.0, 10.0), param)
    solution = solve(problem, RK4(), reltol = 1e-4, abstol = 1e-4)
    data = hcat(solution.u...)
    f₋ = data[(2 + k):(end - 1), end]
    f₊ = data[2:k, end]
    f₀ = data[1, end]
    return (f₀, f₊, f₋)
end
