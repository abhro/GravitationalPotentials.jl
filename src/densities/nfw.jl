# SPDX-License-Identifier: MIT

"""
$(TYPEDEF)

# Fields
$(TYPEDFIELDS)
"""
Base.@kwdef struct NFWDensity{L,MD} <: MassDensityModel
    "Scale density"
    ρ₀::MD
    "Scale radius"
    r₀::L
end

@doc raw"""
    mass_density(model::NFWDensity, s, φ, z)

Density model
```math
ρ(s, φ, z) = \frac{ρ₀}{\frac{r}{r_0} \left(1 + \frac{r}{r_0}\right)^2}
```
where ``r = \sqrt{s^2 + z^2}``
"""
function mass_density(model::NFWDensity{L}, s::L, φ, z::L) where {L}
    r = hypot(s, z)
    return model.ρ₀ / ((r/model.r₀) * (1 + r/model.r₀)^2)
end

Extents.extent(model::NFWDensity{L}) where {L} =
    Extents.Extent((zero(L), L(Inf)), (0, 2π), (L(-Inf), L(Inf)))
