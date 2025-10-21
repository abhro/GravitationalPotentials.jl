# SPDX-License-Identifier: MIT

"""
$(TYPEDEF)

Model for a uniform density cylinder centered at the origin.

### Fields
$(TYPEDFIELDS)

See also [the implementation of `mass_density` for this type](@ref mass_density(::UniformCylinderDensity, ::Any, ::Any, ::Any)).
"""
Base.@kwdef struct UniformCylinderDensity{L,MD} <: MassDensityModel
    "radius of the cylinder"
    R::L
    "height of the cylinder"
    H::L
    "density of the cylinder"
    ρ₀::MD
end
UniformCylinderDensity(R, H, ρ₀) = # for handling mismatching types
    UniformCylinderDensity(promote(R, H)..., ρ₀)


@doc raw"""
    mass_density(model::UniformCylinderDensity, s, φ, z)

Density model
```math
ρ(s, φ, z) = \begin{cases}
    ρ_0 & \text{if } s ≤ R, |z| ≤ H \\
    0 & \text{otherwise}
\end{cases}
```
where
- ``ρ_0`` = `model.ρ₀`
- ``R`` = `model.R`
- ``H`` = `model.H`
"""
function mass_density(model::UniformCylinderDensity{L,MD}, s::L, φ, z::L) where {L,MD}
    return s ≤ model.R && abs(z) ≤ model.H ? model.ρ₀ : zero(MD)
end

mass(model::UniformCylinderDensity) = 2π*model.R^2 * model.H * model.ρ₀

Extents.extent(model::UniformCylinderDensity{L}) where {L} =
    Extents.Extent(s = (zero(L), model.R), φ = (0, 2π), z = (-model.H, model.H))

"""
    onaxispotential(model::UniformCylinderDensity, z)

Potential for test point at ``(0, φ, z)``.
Technically ``φ`` is not well-defined for ``s=0``, but it's irrelevant here.

Mathematical definition
```math
\\Phi(z) = - Gπ ρ_0 \\left[
    (z+h)A - (z-h)B
    + R^2 \\left(
        \\operatorname{artanh}\\left(\\frac{z+h}{A}\\right)
        - \\operatorname{artanh}\\left(\\frac{z-h}{B}\\right)
    \\right)
    - h \\left((2z+h) \\operatorname{sgn}(z+h) + (2z-h) \\operatorname{sgn}(z-h)\\right)
\\right]
```
where ``A = \\sqrt{R^2 + (z+h)^2}``, ``B = \\sqrt{R^2 + (z-h)^2}``, and ``h = H/2``.
"""
function onaxispotential(model::UniformCylinderDensity{L}, z::L) where {L}
    (; R, H) = model
    h = H/2

    # auxiliary quantities, no idea what to call them
    A = hypot(R, z+h)
    B = hypot(R, z-h)

    s1 = sign(z+h)
    s2 = sign(z-h)

    term1 = (z+h) * A - (z-h) * B
    term2 = R^2 * (atanh((z+h)/A) - atanh((z-h)/B))
    term3 = h * ((2z+h) * s1 + (2z-h) * s2)

    return -G * π * model.ρ₀ * (term1 + term2 - term3)
end
