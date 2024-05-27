# SPDX-License-Identifier: MIT

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
    if s ≤ model.r_c && abs(z) ≤ model.h_c
        return model.ρ_c
    end
    return 0.0
end

mass(model::UniformCylinderDensity) = 2π*model.r_c^2 * model.h_c * model.ρ_c
