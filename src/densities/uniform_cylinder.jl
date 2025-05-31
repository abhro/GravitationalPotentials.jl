# SPDX-License-Identifier: MIT

"""
$(TYPEDEF)

Model for a uniform density cylinder centered at the origin.

# Fields
$(TYPEDFIELDS)

See also [the implementation of `mass_density` for this type](@ref mass_density(::UniformCylinderDensity, ::Any, ::Any, ::Any)).
"""
Base.@kwdef struct UniformCylinderDensity{L,MD} <: MassDensityModel
    "radius of the cylinder"
    r_c::L
    "height of the cylinder"
    h_c::L
    "density of the cylinder"
    ρ_c::MD
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
    return s ≤ model.r_c && abs(z) ≤ model.h_c ? model.ρ_c : 0.0
end

mass(model::UniformCylinderDensity) = 2π*model.r_c^2 * model.h_c * model.ρ_c

Extents.extent(model::UniformCylinderDensity) = Extents.Extent(
    s = (0, model.r_c), φ = (0, 2π), z = (-model.h_c, model.h_c))
