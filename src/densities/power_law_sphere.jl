# SPDX-License-Identifier: MIT

"""
$(TYPEDEF)

Model for a sphere with a power law density centered at the origin.

### Fields
$(TYPEDFIELDS)

See also [the implementation of `mass_density` for this type](@ref mass_density(::PowerLawSphereDensity, ::Any, ::Any, ::Any)).
"""
Base.@kwdef struct PowerLawSphereDensity{L,MD} <: MassDensityModel
    "radius of the sphere"
    rₛ::L
    "scale density of the sphere, density changes by `r₀`^`α` at `radius` = `r₀`"
    ρ₀::MD
    "scale radius of the sphere"
    r₀::L
    "power by which density reduces with respect to radius"
    α::Float64

    function PowerLawSphereDensity(rₛ::L, ρ₀::MD, r₀::L, α) where {L,MD}
        if iszero(α)
            throw(DomainError(α, "α cannot be 0. For α = 0 use UniformSphereDensity."))
        end
        return new{L,MD}(rₛ, ρ₀, r₀, α)
    end
end
function PowerLawSphereDensity(rₛ, ρ₀, r₀, α) # for handling mismatching types
    rₛ, r₀ = promote(rₛ, r₀)
    return PowerLawSphereDensity(rₛ, ρ₀, r₀, α)
end

@doc raw"""
    mass_density(model::PowerLawSphereDensity, s, φ, z)

Density model
```math
ρ(s, φ, z) = \begin{cases}
    ρ_0 \left(\sqrt{s^2 + z^2}/R_0\right)^α & \text{if } \sqrt{s^2 + z^2} ≤ R_s \\
    0 & \text{otherwise}
\end{cases}
```
where
- ``ρ_0`` = `model.ρ₀`
- ``R_0`` = `model.r₀`
- ``α`` = `model.α`
- ``R_s`` = `model.rₛ`
"""
function mass_density(model::PowerLawSphereDensity{L,MD}, s::L, φ, z::L) where {L,MD}
    rad = hypot(s, z)
    return rad ≤ model.rₛ ? model.ρ₀ * (rad/model.r₀)^model.α : zero(MD)
end

mass(model::PowerLawSphereDensity) =
    4π/(α+3) * model.rₛ^3 * (model.rₛ/model.r₀)^model.α

Extents.extent(model::PowerLawSphereDensity{L,MD}) where {L,MD} =
    Extents.Extent(s = (zero(L), model.r₀), φ = (0, 2π), z = (-model.r₀, model.r₀))

function potential(model::PowerLawSphereDensity{L}, s::L, φ, z::L) where {L}
    r = hypot(s, z)
    if r ≤ model.rₛ
        error("Not yet implemented")
    end
    M = mass(model)
    return - G * M / r
end
