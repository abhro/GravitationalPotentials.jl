function potential(model::UniformSphereDensity, s, φ, z)
    r = hypot(s, z)
    M = mass(model)

    # outside the sphere, treat it as a point mass
    if r ≥ model.radius
        return -G * M / r
    end
    return -G * M / (2*model.radius) * (3 - r^2/model.radius^2)
end
