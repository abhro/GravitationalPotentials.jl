# SPDX-License-Identifier: MIT

import LinearAlgebra

nonnegative(x::Real) = x < 0 && throw(DomainError(x, "Value must be ≥ 0"))
positive(x::Real) = x ≤ 0 && throw(DomainError(x, "Value must be > 0"))

# follow convention of Griffth's Electrodynamics when naming coordinates
struct Polar{T<:Real,A<:Real}
    s::T
    φ::A
    function Polar(s::T, φ::A) where {T<:Real,A<:Real}
        nonnegative(s)
        return new{T,A}(s, mod2pi(φ))
    end
end
struct Cylindrical{T<:Real,A<:Real}
    s::T
    φ::A # azimuthal angle
    z::T
    function Cylindrical(s::T, φ::A, z::T) where {T<:Real,A<:Real}
        nonnegative(s)
        return new{T,A}(s, mod2pi(φ), z)
    end
end
struct Spherical{T<:Real,A<:Real}
    r::T
    θ::A # polar angle
    φ::A # azimuthal angle
    function Cylindrical(r::T, θ::A, φ::T) where {T<:Real,A<:Real}
        nonnegative(r)
        return new{T,A}(r, mod2pi(θ), mod2pi(φ))
    end
end

## Plane-polar coordinates

LinearAlgebra.norm(x::Polar) = x.s
Base.:-(x::Polar) = Polar(x.s, mod2pi(x.φ + π))

function Base.:+(x₁::Polar, x₂::Polar)
    s₁, φ₁, s₂, φ₂ = x₁.s, x₁.φ, x₂.s, x₂.φ

    s = sqrt(s₁^2 + s₂^2 - 2*s₁*s₂*cos(φ₁-φ₂))
    # TODO simplify/expand?
    φ = atan(s₁*sin(φ₁) + s₂*sin(φ₂), s₁*cos(φ₁) + s₂*cos(φ₂))

    return Polar(s, φ)
end

function Base.:-(x₁::Polar, x₂::Polar)
    s₁, φ₁, s₂, φ₂ = x₁.s, x₁.φ, x₂.s, x₂.φ

    s = sqrt(s₁^2 + s₂^2 + 2*s₁*s₂*cos(φ₁-φ₂))
    # TODO simplify/expand?
    φdiff = atan(s₁*sin(φ₁) - s₂*sin(φ₂), s₁*cos(φ₁) - s₂*cos(φ₂))

    return Polar(s, φdiff)
end

Base.:*(c::Real, x::Polar) = Polar(abs(c) * x.s, x.φ + π*(c<0))
Base.:*(x::Polar, c::Real) = Polar(x.s * abs(c), x.φ + π*(c<0))

Base.:/(x::Polar, c::Real) = Polar(x.s / abs(c), x.φ + π*(c<0))


## Cylindrical coordinates

LinearAlgebra.norm(x::Cylindrical) = hypot(x.s, x.z)
Base.:-(x::Cylindrical) = Cylindrical(x.s, mod2pi(x.φ + π), -x.z)

function Base.:+(x₁::Cylindrical, x₂::Cylindrical)
    s₁, φ₁, z₁, s₂, φ₂, z₂ = x₁.s, x₁.φ, x₁.z, x₂.s, x₂.φ, x₂.z

    ssum = sqrt(s₁^2 + s₂^2 - 2*s₁*s₂*cos(φ₁-φ₂))
    # TODO simplify/expand?
    φsum = atan(s₁*sin(φ₁) + s₂*sin(φ₂), s₁*cos(φ₁) + s₂*cos(φ₂))
    zsum = z₁ + z₂

    return Cylindrical(ssum, φsum, zsum)
end

function Base.:-(x₁::Cylindrical, x₂::Cylindrical)
    s₁, φ₁, z₁, s₂, φ₂, z₂ = x₁.s, x₁.φ, x₁.z, x₂.s, x₂.φ, x₂.z

    sdiff = sqrt(s₁^2 + s₂^2 + 2*s₁*s₂*cos(φ₁-φ₂))
    # TODO simplify/expand?
    φdiff = atan(s₁*sin(φ₁) - s₂*sin(φ₂), s₁*cos(φ₁) - s₂*cos(φ₂))
    zdiff = z₁ - z₂

    return Cylindrical(sdiff, φdiff, zdiff)
end

# the +π trick is to take care of negative c
Base.:*(c::Real, x::Cylindrical) = Cylindrical(abs(c) * x.s, x.φ + π*(c<0), c * x.z)
Base.:*(x::Cylindrical, c::Real) = Cylindrical(x.s * abs(c), x.φ + π*(c<0), x.z * c)

Base.:/(x::Cylindrical, c::Real) = Cylindrical(x.s / abs(c), x.φ + π*(c<0), x.z / c)


## Spherical coordinates
LinearAlgebra.norm(x::Spherical) = x.r
