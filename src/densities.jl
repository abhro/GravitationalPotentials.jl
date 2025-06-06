# SPDX-License-Identifier: MIT

using DocStringExtensions

"""
$(TYPEDEF)

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
function mass_density(model::MassDensityModel, s, φ, z)
    s, z = promote(s, z)
    return mass_density(model, s, φ, z)
end

mass_density(model::MassDensityModel, rvec::NTuple{3,<:Real}) =
    mass_density(model, rvec[1], rvec[2], rvec[3])
mass_density(model::MassDensityModel, rvec::AbstractVector) =
    mass_density(model, rvec...)


include("densities/uniform_sphere.jl")
include("densities/power_law_sphere.jl")
include("densities/uniform_cylinder.jl")
include("densities/power_law_cylinder.jl")
include("densities/spiral_galaxy.jl")
include("densities/einasto.jl")
include("densities/nfw.jl")
