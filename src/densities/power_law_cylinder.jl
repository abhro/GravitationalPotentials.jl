# SPDX-License-Identifier: MIT

"""
$(TYPEDEF)

Model for a cylinder with a power law density centered at the origin.

### Fields
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
function mass_density(model::PowerLawCylinderDensity{L,MD}, s::L, φ, z::L) where {L,MD}
    if s ≤ model.r_c && z ≤ model.h_c
        return model.ρ₀ * (s / model.r₀) ^ model.α * (z / model.h₀)^model.β
    end
    return zero(MD)
end

function mass(model::PowerLawCylinderDensity)
    (; r_c, h_c, ρ₀, r₀, h₉, α, β) = model
    return π * ρ₀ *
           r_c^(α+2) / ((α+2) * r₀^α) *
           h_c^(β+1) / ((β+1) * h₀^β)
end

Extents.extent(model::PowerLawCylinderDensity{L}) where {L} =
    Extents.Extent(s = (zero(L), model.r₀), φ = (0, 2π), z = (-model.h₀, model.h₀))
