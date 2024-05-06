# SPDX-License-Identifier: MIT

abstract type MassDensityModel end

struct UniformSphereDensity <: MassDensityModel
    radius::Float64
    density::Float64
end

struct PowerLawSphereDensity <: MassDensityModel
    radius::Float64
    scale_density::Float64
    scale_radius::Float64
    alpha::Float64
end

struct UniformCylinderDensity <: MassDensityModel
    radius::Float64
    height::Float64
    density::Float64
end

struct SpiralGalaxyDensity <: MassDensityModel
    bulge_radius::Float64
    disk_radius::Float64
    disk_height::Float64

    bulge_density::Float64
    disk_density::Float64
end


function mass_density(model::UniformSphereDensity, s, φ, z)
    rad = hypot(s, z)
    if rad ≤ model.radius
        return model.density
    end
    return 0.0
end

function mass_density(model::PowerLawSphereDensity, s, φ, z)
    rad = hypot(s, z)
    if rad ≤ model.radius
        return model.scale_density * (rad / model.scale_radius) ^ model.alpha
    end
    return 0.0
end

function mass_density(model::UniformCylinderDensity, s, φ, z)
    if s ≤ model.radius && z ≤ model.height
        return model.density
    end
    return 0.0
end

function mass_density(model::SpiralGalaxyDensity, s, φ, z)
    rad = hypot(s, z)

    if rad ≤ model.bulge_radius
        return model.bulge_density
    end
    if s ≤ model.disk_radius && abs(z) ≤ model.disk_height
        return model.disk_density
    end

    return 0.0
end

function mass_density end