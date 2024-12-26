# SPDX-License-Identifier: MIT

module GravitationalPotentials

include("densities.jl")
include("potentials.jl")

export mass_density
export MassDensityModel
export UniformSphereDensity
export PowerLawSphereDensity
export UniformCylinderDensity
export PowerLawCylinderDensity
export SpiralGalaxyDensity
export EinastoDensity


"""Cylindrical coordinates bounds `[(s_lo, s_hi), (φ_lo, φ_hi), (z_lo, z_hi)]`"""
function bounds end
export bounds

"""Cylindrical coordinates bounds `([s_lo, φ_lo, z_lo], [s_hi, φ_hi, z_hi])`"""
function bounds_t(model::MassDensityModel)
    r_bounds, φ_bounds, z_bounds = bounds(model)
    return ([r_bounds[1], φ_bounds[1], z_bounds[1]], [r_bounds[2], φ_bounds[2], z_bounds[2]])
end
export bounds_t


export mass
export potential
end
