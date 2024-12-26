# SPDX-License-Identifier: MIT

"""
    UniformSphereDensity <: MassDensityModel

Model for a uniform density sphere centered at the origin.

# Fields
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
function mass_density(model::UniformSphereDensity, s, φ, z)
    rad = hypot(s, z)
    if rad ≤ model.rₛ
        return model.ρₛ
    end
    return 0.0
end

mass(model::UniformSphereDensity) = 4π/3 * model.ρₛ * model.rₛ^3

function bounds(model::UniformSphereDensity)
    return [(0, model.rₛ), (0, 2π), (-model.rₛ, model.rₛ)]
end
