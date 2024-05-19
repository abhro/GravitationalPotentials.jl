# SPDX-License-Identifier: MIT

"""
    MassDensityModel

Abstract type for all mass density models. Any subtype of `MassDensityModel`
should have a corresponding [`mass_density`](@ref) implementation.
"""
abstract type MassDensityModel end

# the models should be treated as scalars when broadcasting
Base.broadcastable(m::MassDensityModel) = Ref(m)


"""
    mass_density(model::MassDensityModel, s, φ, z)

Given a mass density profile/model, return the corresponding mass density at
cylindrical coordinates ``(s, φ, z)``.
"""
function mass_density end

mass_density(model::MassDensityModel, rvec::NTuple{3,<:Real}) =
    mass_density(model, rvec[1], rvec[2], rvec[3])

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

"""
    PowerLawSphereDensity <: MassDensityModel

Model for a sphere with a power law density centered at the origin.

# Fields
- `rₛ`: radius of the sphere
- `ρ₀`: scale density of the sphere, density changes by `r₀`^`α` at `radius` = `r₀`
- `r₀`: scale radius of the sphere
- `α`: power by which density reduces with respect to radius

See also [the implementation of `mass_density` for this type](@ref mass_density(::PowerLawSphereDensity, ::Any, ::Any, ::Any)).
"""
Base.@kwdef struct PowerLawSphereDensity <: MassDensityModel
    rₛ::Float64
    ρ₀::Float64
    r₀::Float64
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

"""
    UniformCylinderDensity <: MassDensityModel

Model for a uniform density cylinder centered at the origin.

# Fields
- `r_c`: radius of the cylinder
- `h_c`: height of the cylinder
- `ρ_c`: density of the cylinder

See also [the implementation of `mass_density` for this type](@ref mass_density(::UniformCylinderDensity, ::Any, ::Any, ::Any)).
"""
Base.@kwdef struct UniformCylinderDensity <: MassDensityModel
    r_c::Float64
    h_c::Float64
    ρ_c::Float64
end

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
function mass_density(model::UniformCylinderDensity, s, φ, z)
    if s ≤ model.r_h && abs(z) ≤ model.h_c
        return model.ρ_c
    end
    return 0.0
end

"""
    SpiralGalaxyDensity <: MassDensityModel

Model for a spiral galaxy that has both a disk and a central bulge,
with the bulge density dominating in the overlap.

The bulge is modeled as a uniform density sphere,
and the disk is modeled as a uniform density cylinder.

# Fields
- `r_bulge`: radius of the bulge
- `r_disk`: radius of the disk
- `h_disk`: height of the disk
- `ρ_bulge`: density of the bulge
- `ρ_disk`: density of the cylinder

See also [the implementation of `mass_density` for this type](@ref mass_density(::SpiralGalaxyDensity, ::Any, ::Any, ::Any)).
"""
Base.@kwdef struct SpiralGalaxyDensity <: MassDensityModel
    r_bulge::Float64
    r_disk::Float64
    h_disk::Float64

    ρ_bulge::Float64
    ρ_disk::Float64
end

@doc raw"""
    mass_density(model::SpiralGalaxyDensity, s, φ, z)

Density model
```math
ρ(s, φ, z) = \begin{cases}
    ρ_b & \text{if } \sqrt{s^2 + z^2} ≤ R_b \\
    ρ_d & \text{if } s ≤ R_d, |z| ≤ H_d \\
    0 & \text{otherwise}
\end{cases}
```
where
- ``ρ_b`` = `model.ρ_bulge`
- ``ρ_d`` = `model.ρ_disk`
- ``R_b`` = `model.r_bulge`
- ``R_d`` = `model.r_disk`
- ``H_d`` = `model.h_disk`
"""
function mass_density(model::SpiralGalaxyDensity, s, φ, z)
    rad = hypot(s, z)

    if rad ≤ model.r_bulge
        return model.ρ_bulge
    end
    if s ≤ model.disk_radius && abs(z) ≤ model.disk_height
        return model.ρ_disk
    end

    return 0.0
end

"""
    EinastoDensity <: MassDensityModel

# Fields
- `ρ₀`: central density
- `a₀`: harmonic mean radius
- `N`: structural parameter
- `k`: normalizing constant

See [Einasto, J. and Haud, U., “Galactic models with massive corona. I. Method”,
_Astronomy and Astrophysics_, vol. 223, no. 1, pp. 89–94,
1989](https://ui.adsabs.harvard.edu/abs/1989A%26A...223...89E).

See also [the implementation of `mass_density` for this type](@ref
mass_density(::EinastoDensity, ::Any, ::Any, ::Any)).
"""
Base.@kwdef struct EinastoDensity <: MassDensityModel
    ρ₀::Float64
    a₀::Float64
    N::Float64
    k::Float64
end

@doc raw"""
    mass_density(model::EinastoDensity, s, φ, z)

Density model
```math
ρ(s, φ, z) = ρ₀ \exp\left[-(a / k a₀)^{1/N}\right]
```
where
- ``a = √{s² + z²}``
- ``ρ₀`` = `model.ρ₀`
- ``a₀`` = `model.a₀`
- ``k`` = `model.k`
- ``N`` = `model.N`
"""
function mass_density(model::EinastoDensity, s, φ, z)
    a = hypot(s, z)
    return model.ρ₀ * exp(-(a / (model.k * model.a₀)^(1/model.N)))
end
