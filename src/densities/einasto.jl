# SPDX-License-Identifier: MIT

"""
$(TYPEDEF)

# Fields
$(TYPEDFIELDS)

See [Einasto, J. and Haud, U., “Galactic models with massive corona. I. Method”,
_Astronomy and Astrophysics_, vol. 223, no. 1, pp. 89–94,
1989](https://ui.adsabs.harvard.edu/abs/1989A%26A...223...89E).

See also [the implementation of `mass_density` for this type](@ref
mass_density(::EinastoDensity, ::Any, ::Any, ::Any)).
"""
Base.@kwdef struct EinastoDensity{L,MD} <: MassDensityModel
    "central density"
    ρ₀::MD
    "harmonic mean radius"
    a₀::L
    "structural parameter"
    N::Float64
    "normalizing constant"
    k::Float64
end

@doc raw"""
    mass_density(model::EinastoDensity, s, φ, z)

Density model
```math
ρ(s, φ, z) = ρ₀ \exp\left[-\left(\frac{a}{k a₀}\right)^{1/N}\right]
```
where
- ``a = \sqrt{s² + z²}``
- ``ρ₀`` = `model.ρ₀`
- ``a₀`` = `model.a₀`
- ``k`` = `model.k`
- ``N`` = `model.N`
"""
function mass_density(model::EinastoDensity, s, φ, z)
    a = hypot(s, z)
    return model.ρ₀ * exp(-(a / (model.k * model.a₀)^(1/model.N)))
end

Extents.extent(model::EinastoDensity) = Extents.Extent(s = (0, Inf), φ = (0, 2π), z = (-Inf, Inf))
