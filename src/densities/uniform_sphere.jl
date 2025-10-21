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
    R::L
    "density of the sphere"
    ρ₀::MD
end

@doc raw"""
    mass_density(model::UniformSphereDensity, s, φ, z)

Density model
```math
ρ(s, φ, z) = \begin{cases}
    ρ_0 & \text{if } \sqrt{s^2 + z^2} ≤ R \\
    0 & \text{otherwise}
\end{cases}
```
where
- ``ρ_0`` = `model.ρ₀`
- ``R`` = `model.R`
"""
function mass_density(model::UniformSphereDensity{L,MD}, s::L, φ, z::L) where {L,MD}
    r = hypot(s, z)
    return r ≤ model.R ? model.ρ₀ : zero(MD)
end

mass(model::UniformSphereDensity) = 4π/3 * model.ρ₀ * model.R^3

Extents.extent(model::UniformSphereDensity{L,MD}) where {L,MD} =
    Extents.Extent(s = (zero(L), model.R), φ = (0, 2π), z = (-model.R, model.R))

function potential(model::UniformSphereDensity{L}, s::L, φ, z::L) where {L}
    r = hypot(s, z)
    M = mass(model)
    R = model.R

    if r ≥ R
        radial_scale = 1/r # outside the sphere, treat it as a point mass
    else
        radial_scale = (3 - (r/R)^2) / 2R
    end
    return -G * M * radial_scale
end
