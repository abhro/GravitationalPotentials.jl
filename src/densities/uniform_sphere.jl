# SPDX-License-Identifier: MIT

"""
$(TYPEDEF)

Model for a uniform density sphere centered at the origin.

### Fields
$(TYPEDFIELDS)

See also [the implementation of `mass_density` for this type](@ref mass_density(::UniformSphereDensity, ::Any, ::Any, ::Any)).
"""
Base.@kwdef struct UniformSphereDensity{L,MD} <: MassDensityModel
    "radius of the sphere"
    rₛ::L
    "density of the sphere"
    ρₛ::MD
end

@doc raw"""
    mass_density(model::UniformSphereDensity, s, φ, z)

Density model
```math
ρ(s, φ, z) = \begin{cases}
    ρ_s & \text{if } \sqrt{s^2 + z^2} ≤ R_s \\
    0 & \text{otherwise}
\end{cases}
```
where
- ``ρ_s`` = `model.ρₛ`
- ``R_s`` = `model.rₛ`
"""
function mass_density(model::UniformSphereDensity{L,MD}, s::L, φ, z::L) where {L,MD}
    rad = hypot(s, z)
    return rad ≤ model.rₛ ? model.ρₛ : zero(MD)
end

mass(model::UniformSphereDensity) = 4π/3 * model.ρₛ * model.rₛ^3

Extents.extent(model::UniformSphereDensity{L,MD}) where {L,MD} =
    Extents.Extent(s = (zero(L), model.rₛ), φ = (0, 2π), z = (-model.rₛ, model.rₛ))

function potential(model::UniformSphereDensity{L}, s::L, φ, z::L) where {L}
    r = hypot(s, z)
    M = mass(model)
    R = model.rₛ

    # outside the sphere, treat it as a point mass
    radial_scale = r ≥ R ? 1/r : (3/2R - r^2/(2R^3))
    return -G * M * radial_scale
end
