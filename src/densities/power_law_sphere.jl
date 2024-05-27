# SPDX-License-Identifier: MIT

"""
    PowerLawSphereDensity <: MassDensityModel

Model for a sphere with a power law density centered at the origin.

# Fields
$(TYPEDFIELDS)

See also [the implementation of `mass_density` for this type](@ref mass_density(::PowerLawSphereDensity, ::Any, ::Any, ::Any)).
"""
Base.@kwdef struct PowerLawSphereDensity <: MassDensityModel
    "radius of the sphere"
    rₛ::Float64
    "scale density of the sphere, density changes by `r₀`^`α` at `radius` = `r₀`"
    ρ₀::Float64
    "scale radius of the sphere"
    r₀::Float64
    "power by which density reduces with respect to radius"
    α::Float64

    function PowerLawSphereDensity(rₛ, ρ₀, r₀, α)
        if iszero(α)
            throw(DomainError("α cannot be 0. For α = 0 use UniformSphereDensity."))
        end
        return new(rₛ, ρ₀, r₀, α)
    end
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
function mass_density(model::PowerLawSphereDensity, s, φ, z)
    rad = hypot(s, z)
    if rad ≤ model.rₛ
        return model.ρ₀ * (rad / model.r₀) ^ model.α
    end
    return 0.0
end

function mass(model::PowerLawSphereDensity)
    return 4π/(α+3) * (model.radius/model.scale_radius)^model.α * model.radius^3
end
