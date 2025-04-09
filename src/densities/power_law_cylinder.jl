# SPDX-License-Identifier: MIT

"""
$(TYPEDEF)

Model for a cylinder with a power law density centered at the origin.

# Fields
$(TYPEDFIELDS)

See also [the implementation of `mass_density` for this type](@ref mass_density(::PowerLawCylinderDensity, ::Any, ::Any, ::Any)).
"""
Base.@kwdef struct PowerLawCylinderDensity{L,MD} <: MassDensityModel
    "radius of the cylinder"
    r_c::L
    "height of the cylinder"
    h_c::L
    "scale density of the sphere"
    ρ₀::MD
    "scale radius of the cylinder"
    r₀::L
    "scale height of the cylinder"
    h₀::L
    "power by which density changes with respect to radius"
    α::Float64
    "power by which density changes with respect to height"
    β::Float64
end

@doc raw"""
    mass_density(model::PowerLawCylinderDensity, s, φ, z)

Density model
```math
ρ(s, φ, z) = \begin{cases}
    ρ_0 \left(s/R_0\right)^α \left(z/h_0\right)^β & \text{if } s ≤ R_c, \, z ≤ H_c \\
    0 & \text{otherwise}
\end{cases}
```
where
- ``ρ_0`` = `model.ρ₀`
- ``R_0`` = `model.r₀`
- ``H_0`` = `model.h₀`
- ``α`` = `model.α`
- ``β`` = `model.β`
- ``R_c`` = `model.r_c`
- ``H_c`` = `model.h_c`
"""
function mass_density(model::PowerLawCylinderDensity, s, φ, z)
    if s ≤ model.r_c && z ≤ model.h_c
        return model.ρ₀ * (s / model.r₀) ^ model.α * (z / model.h₀)^model.β
    end
    return 0.0
end

function mass(model::PowerLawCylinderDensity)
    α = model.α
    β = model.β
    return π *
           model.r_c^(α+2) / ((α + 2) * model.r₀^α) *
           model.h_c^(β+1) / ((β+1) * mode.h₀^β) *
           model.ρ_0
end

function bounds(model::PowerLawCylinderDensity)
    return [(0, model.r₀), (0, 2π), (-model.h₀, model.h₀)]
end
