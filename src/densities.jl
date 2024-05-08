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
- `radius`: radius of the sphere
- `density`: density of the sphere

See also [the implementation of `mass_density` for this type](@ref mass_density(::UniformSphereDensity, ::Any, ::Any, ::Any)).
"""
Base.@kwdef struct UniformSphereDensity <: MassDensityModel
    radius::Float64
    density::Float64
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
- ``ρ_s`` = `model.density`
- ``R_s`` = `model.radius`
"""
function mass_density(model::UniformSphereDensity, s, φ, z)
    rad = hypot(s, z)
    if rad ≤ model.radius
        return model.density
    end
    return 0.0
end

"""
    PowerLawSphereDensity <: MassDensityModel

Model for a sphere with a power law density centered at the origin.

# Fields
- `radius`: radius of the sphere
- `scale_density`: scale density of the sphere, density changes by
   `scale_density`^`alpha` at `radius` = `scale_radius`
- `scale_radius`: scale radius of the sphere
- `alpha`: power by which density reduces with respect to radius

See also [the implementation of `mass_density` for this type](@ref mass_density(::PowerLawSphereDensity, ::Any, ::Any, ::Any)).
"""
Base.@kwdef struct PowerLawSphereDensity <: MassDensityModel
    radius::Float64
    scale_density::Float64
    scale_radius::Float64
    alpha::Float64
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
- ``ρ_0`` = `model.scale_density`
- ``R_0`` = `model.scale_radius`
- ``α`` = `model.alpha`
- ``R_s`` = `model.radius`
"""
function mass_density(model::PowerLawSphereDensity, s, φ, z)
    rad = hypot(s, z)
    if rad ≤ model.radius
        return model.scale_density * (rad / model.scale_radius) ^ model.alpha
    end
    return 0.0
end

"""
    UniformCylinderDensity <: MassDensityModel

Model for a uniform density cylinder centered at the origin.

# Fields
- `radius`: radius of the cylinder
- `height`: height of the cylinder
- `density`: density of the cylinder

See also [the implementation of `mass_density` for this type](@ref mass_density(::UniformCylinderDensity, ::Any, ::Any, ::Any)).
"""
Base.@kwdef struct UniformCylinderDensity <: MassDensityModel
    radius::Float64
    height::Float64
    density::Float64
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
- ``ρ_0`` = `model.density`
- ``R_c`` = `model.radius`
- ``H_c`` = `model.height`
"""
function mass_density(model::UniformCylinderDensity, s, φ, z)
    if s ≤ model.radius && abs(z) ≤ model.height
        return model.density
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
- `bulge_radius`: radius of the bulge
- `disk_radius`: radius of the disk
- `disk_height`: height of the disk
- `bulge_density`: density of the bulge
- `disk_density`: density of the cylinder

See also [the implementation of `mass_density` for this type](@ref mass_density(::SpiralGalaxyDensity, ::Any, ::Any, ::Any)).
"""
Base.@kwdef struct SpiralGalaxyDensity <: MassDensityModel
    bulge_radius::Float64
    disk_radius::Float64
    disk_height::Float64

    bulge_density::Float64
    disk_density::Float64
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
- ``ρ_b`` = `model.bulge_density`
- ``ρ_d`` = `model.disk_density`
- ``R_b`` = `model.bulge_radius`
- ``R_d`` = `model.disk_radius`
- ``H_d`` = `model.disk_height`
"""
function mass_density(model::SpiralGalaxyDensity, s, φ, z)
    rad = hypot(s, z)

    if rad ≤ model.bulge_radius
        return model.bulge_density
    end
    if s ≤ model.disk_radius && abs(z) ≤ model.disk_height
        return model.disk_density
    end

    return 0.0
end