abstract type CarrierGas end
struct Air <: CarrierGas end
struct N2 <: CarrierGas end

const kb = 1.38064852e-23    # Boltzmann constant      [J K-1 molec-1]
const r = 8.314472           # Universal gas const     [J K-1 mol-1]
const ec = 1.602176565e-19   # Elemental charge        [C]
const eps = 8.854187817e-12  # Dielectric constant     [farad m-1]
const na = 6.02214086e23     # Avagadro's constant     [molecule mol-1]
const mwair = 28.9647e-3     # Molecular weight of air [kg mol-1]

@doc raw"""
    DMAconfig

Data type to abstract the DMA geometry and state of the fluid.

    t::AbstractFloat          # Temperature [K]
    p::AbstractFloat          # Pressure [Pa]
    qsa::AbstractFloat        # Sample flow [m3 s-1]
    qsh::AbstractFloat        # Sheath flow [m3 s-1]
    r1::AbstractFloat         # Inner column radius [m]
    r2::AbstractFloat         # Outer column radius [m]
    l::AbstractFloat          # Column length [m]
    leff::AbstractFloat       # Effective length [m]
    polarity::Symbol          # Power supply polarity [:+] or [:-]
    m::Int8                   # Number of charges in charge correction [-]
    DMAtype::Symbol           # Designate :radial, :cylindrical

Example Usage
```julia
t,p = 295.15, 1e5
qsa,qsh = 1.66e-5, 8.3e-5
r₁,r₂,l = 9.37e-3,1.961e-2,0.44369
Λ = DMAconfig(t,p,qsa,qsh,r₁,r₂,l,13.0,:-,6,:cylindrical)
```julia

!!! note

    When defining a radial DMA, r₁,r₂,l map to  r₁,r₂,b as defined in Zhang Shou-Hua Zhang,
    Yoshiaki Akutsu, Lynn M. Russell, Richard C. Flagan & John H. Seinfeld (1995)
    Radial Differential Mobility Analyzer, Aerosol Science and Technology,
    23:3, 357-372, DOI: 10.1080/02786829508965320.

"""
struct DMAconfig
    t::AbstractFloat          # Temperature [K]
    p::AbstractFloat          # Pressure [Pa]
    qsa::AbstractFloat        # Sample flow [m3 s-1]
    qsh::AbstractFloat        # Sheath flow [m3 s-1]
    r1::AbstractFloat         # Inner column radius [m]
    r2::AbstractFloat         # Outer column radius [m]
    l::AbstractFloat          # Column length [m]
    leff::AbstractFloat       # Effective length [m]
    polarity::Symbol          # Power supply polarity [:+] or [:-]
    m::Int8                   # Number of charges in charge correction [-]
    DMAtype::Symbol           # Designate :radial, :cylindrical
    gas::CarrierGas           # CarrierGas, either Air() or N2()
end

DMAconfig(t, p, qsa, qsh, r1, r2, l, leff, polarity, m, DMAtype) =
    DMAconfig(t, p, qsa, qsh, r1, r2, l, leff, polarity, m, DMAtype, Air())

"""
    clean(x)

Defined as shorthand:

```julia
clean(x) = map(x -> x < 0.0 ? 0.0 : x, x)
```

The function removes negative numbers and set them zero. It is used to cleanup
inverted size distribution data, which may contain small negative values from
inversion noise.
"""
clean(x) = map(x -> x < 0.0 ? 0.0 : x, x)

@doc raw"""
    λ(Λ:DMAconfig, ::CarrierGas)

λ is the mean free path of air in [m]. It  depends on temperature [K] and Pressure [Pa].
Temperature and pressure are taken from the DMA configuration. Currently dry air
and N2 are supported.

``\lambda = 6.6 \times 10^{-8}\frac{101315}{p}\frac{t}{293.15}``

Example Usage

```julia
t,p = 295.15, 1e5
qsa,qsh = 1.66e-5, 8.3e-5
r₁,r₂,l = 9.37e-3,1.961e-2,0.44369
Λ = DMAconfig(t,p,qsa,qsh,r₁,r₂,l,0.0,:-,6,:cylindrical,N2())
mfp = λ(Λ, N2())
```
"""
λ(Λ::DMAconfig, ::Air) = 6.6e-8 * (101315.0 ./ Λ.p) * (Λ.t / 293.15)
λ(Λ::DMAconfig, ::N2) = 5.9e-8 * (101315.0 ./ Λ.p) * (Λ.t / 293.15)

@doc raw"""
    η(Λ::DMAconfig, ::CarrierGas)

η is the viscosity of air in [Pa s] and depends on temperature [K]. Temperature
is taken from the DMA configuration. Currently dry air and N2 is supported.

``\eta = 1.83245\times10^{-5} \exp \left(1.5 \ln \left[\frac{T}{296.1}\right]\right)\left
(\frac{406.55}{T+110.4} \right)``

Example Usage

```julia
t,p = 295.15, 1e5
qsa,qsh = 1.66e-5, 8.3e-5
r₁,r₂,l = 9.37e-3,1.961e-2,0.44369
Λ = DMAconfig(t,p,qsa,qsh,r₁,r₂,l,0.0,:-,6,:cylindrical,Air())
viscosity = η(Λ, Air())
```
"""
η(Λ::DMAconfig, ::Air) = 1.83245e-5 * exp(1.5 * log(Λ.t / 296.1)) * (406.55) / (Λ.t + 110.4)
η(Λ::DMAconfig, ::N2) = 1.663e-5 * (Λ.t / 273.0)^1.5 * (380.0) / (Λ.t + 107.0)

@doc raw"""
    cc(Λ::DMAconfig, d)

Cunningham slip-flow correction factor. The slip flow correction accounts for the
decreased drag particles experience relative to Stokes' drag force when particle
size approaches the scale of the mean free path of air. It is computed following
Hinds (1999) Eq. 3.20. Temperature and pressure are taken from the DMA configuration.
The units of diameter are in [m] and the function accepts scalars or arrays.

``
c_c = 1+\frac{\lambda}{d_p} \left(2.34+1.05 \exp \left[-0.39 \frac{d_p}{\lambda}\right]\right)
``

where ``d_p`` is the particle diameter and ``\lambda`` is the mean
free path of air, which is computed as a function of pressure and temperature.

Example Usage

```julia
Dp = exp10.(range(log10(1e-9), stop=log10(1000e-9), length=100))
t,p = 295.15, 1e5
qsa,qsh = 1.66e-5, 8.3e-5
r₁,r₂,l = 9.37e-3,1.961e-2,0.4436, Λ.gas9
Λ = DMAconfig(t,p,qsa,qsh,r₁,r₂,l,0.0,:-,6,:cylindrical)
correction = cc(Λ, Dp)
```
"""
cc(Λ::DMAconfig, d) =
    1.0 .+ λ(Λ, Λ.gas) ./ d .* (2.34 .+ 1.05 .* exp.(-0.39 .* d ./ λ(Λ, Λ.gas)))

@doc raw"""
    dab(Λ::DMAconfig, d)

The diffusion coefficient of particles in air, ``d_{ab}``, describes the random
displacement of particles due to Brownian motion. It is computed via the Stokes-Einstein
relation (Hinds, 1999, Eq. 7.20). Temperature and pressure are taken from the DMA
configuration. The units of diameter are in [m] and the function accepts scalars or arrays.

``d_{ab} = \frac{k_bTc_c}{3\pi\eta d_p}``

where ``k_b`` is Boltzmann's constant and ``\eta`` is the viscosity of air in [Pa s],
``c_c`` is the Cunningham slip flow correction and ``d_p`` is the particle diameter.
``d_{ab}`` is in [m² s⁻¹].

Example Usage

```julia
Dp = exp10.(range(log10(1e-9), stop=log10(1000e-9), length=100))
t,p = 295.15, 1e5
qsa,qsh = 1.66e-5, 8.3e-5
r₁,r₂,l = 9.37e-3,1.961e-2,0.44369
Λ = DMAconfig(t,p,qsa,qsh,r₁,r₂,l,0.0,:-,6,:cylindrical)
diffusion_coefficient = dab(Λ,Dp)
```
"""
dab(Λ::DMAconfig, d) = kb * Λ.t * cc(Λ, d) ./ (3.0π * η(Λ, Λ.gas) * d)

u = (Λ, D) -> @. D * Λ.leff / Λ.qsa
Peff =
    u -> @. 0.82 * exp(-11.5u) +
       0.1 * exp(-70.0u) +
       0.03 * exp(-180.0u) +
       0.02 * exp(-340.0u)
Taf = (x, μ, σ) -> @. 0.5 * (1.0 + erf((x - μ) / (√2σ)))

@doc raw"""
    Tl(Λ::DMAconfig, Dp)

Penetration efficiency through the TSI cylindrical DMA using the parameterization by
Reineking & Porstendörfer (1986). The particle diameter Dp is in [nm].

``T_l = 0.82\exp(-11.5u)+0.1\exp(-70.0u)+0.03\exp(-180.0u)+0.02\exp(-340.0u)``

where ``u = \frac{d_{ab} l_{eff}}{q_{sa}}$, $l_{eff}`` is the parameterized effective
diffusion length, and ``q_{sa}`` is the aerosol flow rate through the DMA.

!!! note

    Λ contains the effective length, aerosol flow rate, temperature and pressure to
    compute ``d_{ab}``. To treat multiple DMAs with different {leff, qsa,
    t, p} in a single script, the function Tl is  embedded  in the
    DifferentialMobilityAnalyzer data type.

"""
function Tl(Λ::DMAconfig, Dp)
    if Λ.leff == 0.0
        return ones(length(Dp))
    else
        return clean(Peff(u(Λ, dab(Λ, Dp * 1e-9))))
    end
end

function Tl(Λ::DMAconfig, Z, k)
    if Λ.leff == 0.0
        return ones(length(Z))
    else
        Dp = map(zs -> ztod(Λ, k, zs) * 1e-9, Z)
        return clean(Peff(u(Λ, dab(Λ, Dp))))
    end
end


@doc raw"""
    dtoz(Λ::DMAconfig, d)

The function returns the mobility ``z`` according to

``d_p =  \frac{kec_c}{3\pi \eta z^s}``

where ``e`` is the elementary charge,  ``k`` is the number of charges on the particle,
``c_c`` is the Cunningham correction factor, and ``\eta`` is the viscosity of the fluid.
The diameter in dtoz is in units of [m].

Example Usage
```julia
t,p = 295.15, 1e5
qsa,qsh = 1.66e-5, 8.3e-5
r₁,r₂,l = 9.37e-3,1.961e-2,0.44369
Λ = DMAconfig(t,p,qsa,qsh,r₁,r₂,l,0.0,:-,6,:cylindrical)
mobility = dtoz(Λ,dp*1e-9) # [m2 V-1 s-1]
```
"""
dtoz(Λ::DMAconfig, d) = dtoz(Λ, 1, d)

dtoz(Λ::DMAconfig, k, d) = k .* ec .* cc(Λ, d) ./ (3.0π .* η(Λ, Λ.gas) .* d)

@doc raw"""
    vtoz(Λ::DMAconfig, v)

Converts between voltage and selected mobility. It is the inverse of [ztov](@ref).

For the cylindrical DMA and balanced flows:

``z^s = \frac{q_{sh}}{2\pi l v} \ln \left(\frac{r_2}{r_1}\right)``

For the radial DMA and balanced flows:

``z^s = \frac{q_{sh} l}{\pi v \left({r_2}^2 - {r_1}^2\right)} ``

where ``v`` is the potential applied between the inner and out section of the annulus,
``r_1``, ``r_2``, and ``l`` are the dimensions of the cylindrical DMA  and ``q_{sh}`` is
the sheath flow rate.

Example Usage
```julia
t,p = 295.15, 1e5
qsa,qsh = 1.66e-5, 8.3e-5
r₁,r₂,l = 9.37e-3,1.961e-2,0.44369
Λ = DMAconfig(t,p,qsa,qsh,r₁,r₂,l,0.0,:-,6,:cylindrical)
mobility = vtoz(Λ,1000.0) # [m2 V-1 s-1]
```
"""
vtoz(Λ::DMAconfig, v) =
    (Λ.DMAtype == :radial) ? Λ.qsh .* Λ.l / (π .* (Λ.r2^2.0 - Λ.r1^2) .* v) :
    Λ.qsh ./ (2.0π .* Λ.l .* v) .* log(Λ.r2 / Λ.r1)

@doc raw"""
    ztov(Λ::DMAconfig, v)

Converts between selected mobility and voltage. It is the inverse of [vtoz](@ref).

Example Usage
```julia
t,p = 295.15, 1e5
qsa,qsh = 1.66e-5, 8.3e-5
r₁,r₂,l = 9.37e-3,1.961e-2,0.44369
Λ = DMAconfig(t,p,qsa,qsh,r₁,r₂,l,0.0,:-,6,:cylindrical)
voltage = ztov(Λ,1e-9)
```
"""
const di = 1e-7
ztov(Λ::DMAconfig, z) =
    (Λ.DMAtype == :radial) ? Λ.qsh .* Λ.l / (π .* (Λ.r2^2.0 - Λ.r1^2) * z) :
    Λ.qsh ./ (2.0π .* Λ.l .* z) .* log(Λ.r2 / Λ.r1)

f(Λ, i, z, di) = @. i .* ec .* cc($Ref(Λ), di) ./ (3.0π .* η($Ref(Λ), $Ref(Λ.gas)) .* z)
converge(f, g) = maximum(abs.(1.0 .- f ./ g) .^ 2.0) < 1e-24
g(Λ, i, z, di) = converge(f(Λ, i, z, di), di) ? di : g(Λ, i, z, f(Λ, i, z, di))

@doc raw"""
    ztod(Λ::DMAconfig, i::Int, z)

Converts mobility to diameter.
- Λ is the DMA configuration
- i is the number of charges
- z is the mobility

```julia
t,p = 295.15, 1e5
qsa,qsh = 1.66e-5, 8.3e-5
r₁,r₂,l = 9.37e-3,1.961e-2,0.44369
Λ = DMAconfig(t,p,qsa,qsh,r₁,r₂,l,0.0,:-,6,:cylindrical)
z = dtoz(Λ,100.0*1e-9)
diameter = ztod(Λ,1,z)
```
"""
@memoize ztod(Λ::DMAconfig, i::Int, z) = g(Λ, i, z, di) .* 1e9;

@doc raw"""
    Ω(Λ::DMAconfig, Z, zs)

The DMA transfer function is the probability that a particle of a particle of a given size
exits the classifier via the sample flow. The diffusive broadened DMA transfer function is
computed assuming blanced sheath and excess flows using the expression of Stolzenburg
and McMurry (2008).

``\Omega(\tilde{z},\beta,\sigma) = \frac{\sigma}{\sqrt{2}\beta}\left[\epsilon \left(
    \frac{\tilde{z}-(1+\beta)}{\sqrt{2}\sigma} \right) + \epsilon \left (\frac{\tilde{z}-
    (1-\beta)}{\sqrt{2}\sigma} \right) - 2\epsilon \left
    ( \frac{\tilde{z}-1}{\sqrt{2}\sigma}\right)  \right]``

where ``\tilde{z} = \frac{z}{z^s}`` is the dimensionless mobility, ``z`` is the particle
mobility ``z^s`` is the centroid mobility selected by the DMA,
``\epsilon = x \mathrm{erf}(x) +\left(\exp(-x^2)/\sqrt{\pi}\right)``, ``\mathrm{erf}`` is
the error function, and ``\beta = \frac{q_{sa}}{q_{sh}}``. The parameter ``\sigma``
accounts for diffusional broading of the transfer function. Assuming plug flow,
``\sigma`` can be computed using the following equations Hagwood (1999)

``\gamma = \left(\frac{r_1}{r_2}\right)^2``

``I = \frac{1}{2}(1+γ)``

``\kappa = \frac{lr_2}{r_2^2-r_1^2}``

``G = \frac{4(1+\beta)^2}{(1-γ)} \left[I+\{2(1+\beta)\kappa\}^{-2} \right ]``

``\sigma = \sqrt{\frac{2G\pi ld_{ab}}{q_{sh}}}``

Inputs for flow are taken from the DMAconfig. The function expects a mobility scalar z or vector Z,
and a centroid mobility zˢ.

Example Usage
```julia
zˢ = dtoz(Λ, 200e-9)      # centroid mobility for Dp = 200 nm
z = [1e-9, 1e-8, 1e-7]    # mobility
Ω(Λ,z,zˢ)                 # Output of the transfer function
```

!!! note

    The function Ω is embedded in the the Type DifferentialMobilityAnalyzers.jl, which
    assigns δ.Ω either to this function Ω or Ωav applicable to scanning mode,

"""
function Ω(Λ::DMAconfig, Z, zs)
    ε = (x) -> @. x * erf(x) .+ exp(-x^2.0) / √π
    D = dab(Λ, ztod(Λ, 1, zs) * 1e-9)
    β = Λ.qsa / Λ.qsh
    γ = (Λ.r1 / Λ.r2)^2.0
    I = 0.5(1.0 + γ)
    κ = Λ.l * Λ.r2 / (Λ.r2^2.0 - Λ.r1^2.0)
    G = 4.0(1.0 + β)^2.0 / (1.0 - γ) * (I + (2.0(1.0 + β)κ)^(-2.0))
    σ = √(G * 2.0π * Λ.l * D / Λ.qsh)
    f =
        (Z, σ, β, ε) -> @. σ / (√2.0 * β) * (
            ε((Z - (1.0 + β)) / (√2.0 * σ)) + ε((Z - (1.0 - β)) / (√2.0 * σ)) -
            2.0 * ε((Z - 1.0) / (√2.0 * σ))
        )
    return clean(f(Z / zs, σ, β, ε))
end


function Ωav(Λ::DMAconfig, Z, zs)
    Vex = mylogspace(Ve[i], Ve[i + 1], nint)
    return mapreduce(zˢ -> Ω(Λ, Z, zˢ / k, k), +, vtoz(Λ, Vex)) / nint
end
