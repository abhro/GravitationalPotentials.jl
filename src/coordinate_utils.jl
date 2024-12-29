# SPDX-License-Identifier: MIT

import LinearAlgebra
using CoordinateTransformations

## Plane-polar coordinates

LinearAlgebra.norm(x::Polar) = x.r
Base.:-(x::Polar) = Polar(x.r, mod2pi(x.θ + π))

function Base.:+(x₁::Polar, x₂::Polar)
    r₁, θ₁, r₂, θ₂ = x₁.r, x₁.θ, x₂.r, x₂.θ

    rsum = sqrt(r₁^2 + r₂^2 - 2*r₁*r₂*cos(θ₁-θ₂))
    # TODO simplify/expand?
    θsum = atan(r₁*sin(θ₁) + r₂*sin(θ₂), r₁*cos(θ₁) + r₂*cos(θ₂))

    return Polar(rsum, θsum)
end

function Base.:-(x₁::Polar, x₂::Polar)
    r₁, θ₁, r₂, θ₂ = x₁.r, x₁.θ, x₂.r, x₂.θ

    rdiff = sqrt(r₁^2 + r₂^2 + 2*r₁*r₂*cos(θ₁-θ₂))
    # TODO simplify/expand?
    θdiff = atan(r₁*sin(θ₁) - r₂*sin(θ₂), r₁*cos(θ₁) - r₂*cos(θ₂))

    return Polar(rdiff, θdiff)
end

Base.:*(c::Number, x::Polar) = Polar(c * x.r, x.θ)
Base.:*(x::Polar, c::Number) = Polar(x.r * c, x.θ)

Base.:/(x::Polar, c::Number) = Polar(x.r / c, x.θ)


## Cylindrical coordinates

LinearAlgebra.norm(x::Cylindrical) = hypot(x.r, x.z)
Base.:-(x::Cylindrical) = Cylindrical(x.r, mod2pi(x.θ + π), -x.z)

function Base.:+(x₁::Cylindrical, x₂::Cylindrical)
    r₁, θ₁, z₁, r₂, θ₂, z₂ = x₁.r, x₁.θ, x₁.z, x₂.r, x₂.θ, x₂.z

    rsum = sqrt(r₁^2 + r₂^2 - 2*r₁*r₂*cos(θ₁-θ₂))
    # TODO simplify/expand?
    θsum = atan(r₁*sin(θ₁) + r₂*sin(θ₂), r₁*cos(θ₁) + r₂*cos(θ₂))
    zsum = z₁ + z₂

    return Cylindrical(rsum, θsum, zsum)
end

function Base.:-(x₁::Cylindrical, x₂::Cylindrical)
    r₁, θ₁, z₁, r₂, θ₂, z₂ = x₁.r, x₁.θ, x₁.z, x₂.r, x₂.θ, x₂.z

    rdiff = sqrt(r₁^2 + r₂^2 + 2*r₁*r₂*cos(θ₁-θ₂))
    # TODO simplify/expand?
    θdiff = atan(r₁*sin(θ₁) - r₂*sin(θ₂), r₁*cos(θ₁) - r₂*cos(θ₂))
    zdiff = z₁ - z₂

    return Cylindrical(rdiff, θdiff, zdiff)
end

# the +π trick is to take care of negative c
Base.:*(c::Real, x::Cylindrical) = Cylindrical(
    abs(c) * x.r, mod2pi(x.θ + π*(c<0)), c * x.z)
Base.:*(x::Cylindrical, c::Real) = Cylindrical(
    x.r * abs(c), mod2pi(x.θ + π*(c<0)), x.z * c)

Base.:/(x::Cylindrical, c::Real) = Cylindrical(
    x.r / abs(c), mod2pi(x.θ + π * (c < 0)), x.z / c)


## Spherical coordinates
LinearAlgebra.norm(x::Spherical) = x.r
