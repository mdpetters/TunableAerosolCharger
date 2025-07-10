# ---------------------------------------------------------------------------------
# electricField.jl
#
# Functions to calculate the electric field inside the charger
#
# Author: Markus Petters (markus.petters@ucr.edu)
# ---------------------------------------------------------------------------------

"""
    hod(f, x, n::Int)

Purpose:

Evaluate the nth derivative of function f using forward mode automatic differentiation.

Inputs
    - f: function
    - x: value at which function is evaluated
    - n: derivative order (n = 1 is first derivate, n = 2 is second derivative

Outputs:
    - n'th order derivative of function f at point x

Example:

The third derivative of f(x) is 24x + 6

```julia
julia> f(x) = x^4 + 2.0*x^3 + x + 1
f (generic function with 1 method)

julia> TunableAerosolCharger.hod(f, 0.5, 3)
24.0
```
"""
function hod(f, x, n::Int)
    if n == 1
        return ForwardDiff.derivative(f, x)
    else
        return ForwardDiff.derivative(x -> hod(f, x, n - 1), x)
    end
end

function Ψ(lens::einzelLens, z, r; maxdiff = 4)
    dz = lens.D / 100.0  # discretization for derivative
    zx = [z - 2 * dz z - dz z z + dz z + 2 * dz]  # 5 step derivative

    # Solve Eq. 4 in Rahid et al. (2010) for centerline
    function psi(z)
        zp = @. (2.0 * z + lens.S) / lens.D
        zm = @. (2.0 * z - lens.S) / lens.D
        Ψ = @. (lens.V2 - lens.V1) * lens.D / (2.0 * pi * lens.S) *
           (zp * atan(zp) - zm * atan(zm)) + (lens.V1 + lens.V2) / 2
    end

    # Eq. (2) in Rahid et al.
    Ψ = psi(z)
    out = map(1:maxdiff) do n
        (-1)^n / factorial(n)^2 * (r / 2)^(2 * n) * hod(psi, z, 2 * n)
    end
    return psi(z) + sum(out)
end

function lensFields(lens::einzelLens, z, r)
    func(x, y) = Ψ(lens, x, y)
    func(v) = func(v...)
    ψ = [func(x, y) for x in z, y in r]
    Ez = [ForwardDiff.gradient(func, [x, y])[1] for x in z, y in r]
    Er = [ForwardDiff.gradient(func, [x, y])[2] for x in z, y in r]
    return ψ', Ez', Er', z, r
end
