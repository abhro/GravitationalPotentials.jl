# SPDX-License-Identifier: MIT

module GravitationalPotentials

using Coordinates
import Extents
using PhysicalConstants.CODATA2022: G

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
export NFWDensity

end
