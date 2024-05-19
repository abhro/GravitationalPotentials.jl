# SPDX-License-Identifier: MIT

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
    if s ≤ model.r_disk && abs(z) ≤ model.h_disk
        return model.ρ_disk
    end

    return 0.0
end
