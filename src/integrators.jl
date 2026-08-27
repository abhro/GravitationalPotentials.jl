# SPDX-License-Identifier: MIT

"""
$(SIGNATURES)

Integrate a function ``f(s)`` over a annulus of inner radius `s1` and outer
radius `s2` centered at the origin. Assumes ``f`` has azimuthal symmetry.
"""
function integrate_disk(f, s1, s2)
    I = quadgk(s -> f(s) * s, s1, s2) |> first
    return 2π * I
end
integrate_disk(f, S) = integrate_disk(f, zero(S), S)

"""
$(SIGNATURES)

Integrate a function ``f(s, z)`` over a cylinder. Assumes azimuthal symmetry.
"""
function integrate_cylinder(f, S, z1, z2)
    z_integrand = z -> first(quadgk(s -> s * f(s, z), zero(S), S))
    I, resid = quadgk(z_integrand, z1, z2)
    return 2π * I
end
