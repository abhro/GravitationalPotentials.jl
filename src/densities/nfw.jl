# SPDX-License-Identifier: MIT

"""
    NFWDensity <: MassDensityModel

# Fields
$(TYPEDFIELDS)
"""
struct NFWDensity{L,MD} <: MassDensityModel
    "Scale density"
    ρ₀::MD
    "Scale radius"
    r_s::L
end

@doc raw"""
    mass_density(model::NFWDensity, s, φ, z)

Density model
```math
ρ(s, φ, z) = \frac{ρ₀}{\frac{r}{r_s} \left(1 + \frac{r}{r_s}\right)^2}
```
where ``r = \sqrt{s^2 + z^2}``
"""
function mass_density(model::NFWDensity, s, φ, z)
    r = hypot(s, z)
    return model.ρ₀ / ((r/model.r_s) * (1 + r/model.r_s)^2)
end

function bounds(model::NFWDensity)
    return [(0, Inf), (0, 2π), (-Inf, Inf)]
end
