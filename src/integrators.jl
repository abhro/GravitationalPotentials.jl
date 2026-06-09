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
integrate_disk(f, s) = integrate_disk(f, zero(s), s)

"""
$(SIGNATURES)

Integrate a function ``f(s, z)`` over a cylinder.
"""
function integrate_cylinder(f, s2, z1, z2)
