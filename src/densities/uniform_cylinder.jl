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
    r_c::L
    "height of the cylinder"
    h_c::L
    "density of the cylinder"
    ρ_c::MD
end
UniformCylinderDensity(r_c, h_c, ρ_c) = # for handling mismatching types
    UniformCylinderDensity(promote(r_c, h_c)..., ρ_c)


@doc raw"""
    mass_density(model::UniformCylinderDensity, s, φ, z)

Density model
```math
ρ(s, φ, z) = \begin{cases}
    ρ_0 & \text{if } s ≤ R_c, |z| ≤ H_c \\
    0 & \text{otherwise}
\end{cases}
```
where
- ``ρ_0`` = `model.ρ_c`
- ``R_c`` = `model.r_c`
- ``H_c`` = `model.h_c`
"""
function mass_density(model::UniformCylinderDensity{L,MD}, s::L, φ, z::L) where {L,MD}
    return s ≤ model.r_c && abs(z) ≤ model.h_c ? model.ρ_c : zero(MD)
end

mass(model::UniformCylinderDensity) = 2π*model.r_c^2 * model.h_c * model.ρ_c

Extents.extent(model::UniformCylinderDensity{L}) where {L} =
    Extents.Extent(s = (zero(L), model.r_c), φ = (0, 2π), z = (-model.h_c, model.h_c))

"""
    onaxispotential(model::UniformCylinderDensity, z)

Potential for test point at ``(0, φ, z)``.
Technically ``φ`` is not well-defined for ``s=0``, but it's irrelevant here.
"""
function onaxispotential(model::UniformCylinderDensity{L}, z::L) where {L}
    R = model.r_c
    H = model.h_c/2

    # auxiliary quantities, no idea what to call them
    A = hypot(R, z+H)
    B = hypot(R, z-H)

    s1 = sign(z+H)
    s2 = sign(z-H)

    term1 = (z+H) * A
    term2 = (z-H) * B
    term3 = R^2 * (atanh((z+H)/A) - atanh((z-H)/B))
    term4 = H * ((2z+H) * s1 + (2z-H) * s2)

    return -G * π * model.ρ_c * (term1 - term2 + term3 - term4)
end
