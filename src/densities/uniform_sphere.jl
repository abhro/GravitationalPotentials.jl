# SPDX-License-Identifier: MIT

"""
    UniformSphereDensity <: MassDensityModel

Model for a uniform density sphere centered at the origin.

# Fields
- `rₛ`: radius of the sphere
- `ρₛ`: density of the sphere

See also [the implementation of `mass_density` for this type](@ref mass_density(::UniformSphereDensity, ::Any, ::Any, ::Any)).
"""
Base.@kwdef struct UniformSphereDensity <: MassDensityModel
    rₛ::Float64
    ρₛ::Float64
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
