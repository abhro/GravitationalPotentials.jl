# SPDX-License-Identifier: MIT

"""
$(TYPEDEF)

Model for a spiral galaxy that has both a disk and a central bulge,
with the bulge density dominating in the overlap.

The bulge is modeled as a uniform density sphere,
and the disk is modeled as a uniform density cylinder.

### Fields
$(TYPEDFIELDS)

See also [the implementation of `mass_density` for this type](@ref mass_density(::SpiralGalaxyDensity, ::Any, ::Any, ::Any)).
"""
Base.@kwdef struct SpiralGalaxyDensity{L,MD} <: MassDensityModel
    "radius of the bulge"
    r_bulge::L
    "radius of the disk"
    r_disk::L
    "height of the disk"
    h_disk::L

    "density of the bulge"
    ρ_bulge::MD
    "density of the disk"
    ρ_disk::MD
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
function mass_density(model::SpiralGalaxyDensity{L,MD}, s::L, φ, z::L) where {L,MD}
    rad = hypot(s, z)

    if rad ≤ model.r_bulge # within the bulge
        return model.ρ_bulge
    end
    if s ≤ model.r_disk && abs(z) ≤ model.h_disk # within the disk
        return model.ρ_disk
    end

    return zero(MD) # outside the galaxy
end

function Extents.extent(model::SpiralGalaxyDensity{L}) where {L}
    (; r_bulge, r_disk, h_disk) = model
    return Extents.Extent(
        s = (zero(L), max(r_bulge, r_disk, h_disk)),
        φ = (0, 2π),
        z = (-h_disk, h_disk))
end
