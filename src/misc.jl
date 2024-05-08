function potential(model::UniformSphereDensity, s, φ, z)
    r = hypot(s, z)
    M = mass(model)

    # outside the sphere, treat it as a point mass
    if r ≥ model.radius
        return -G * M / r
    end
    return -G * M / (2*model.radius) * (3 - r^2/model.radius^2)
end

mass(model::UniformSphereDensity) = 4π/3 * model.density * model.radius^3

function mass(model::PowerLawSphereDensity)
    α = model.alpha
    return 4π/(α+3) * (model.radius/model.scale_radius)^α * model.radius^3
end

mass(model::UniformCylinderDensity) = 2π*model.radius^2 * model.height * model.density