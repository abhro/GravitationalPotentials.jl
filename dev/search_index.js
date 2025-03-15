var documenterSearchIndex = {"docs":
[{"location":"api/#API-reference","page":"API reference","title":"API reference","text":"","category":"section"},{"location":"api/","page":"API reference","title":"API reference","text":"Documentation for GravitationalPotentials.","category":"page"},{"location":"api/#Index","page":"API reference","title":"Index","text":"","category":"section"},{"location":"api/","page":"API reference","title":"API reference","text":"","category":"page"},{"location":"api/#Docstrings","page":"API reference","title":"Docstrings","text":"","category":"section"},{"location":"api/","page":"API reference","title":"API reference","text":"bounds\nbounds_t","category":"page"},{"location":"api/#GravitationalPotentials.bounds","page":"API reference","title":"GravitationalPotentials.bounds","text":"bounds(model)\n\nCylindrical coordinates bounds [(s_lo, s_hi), (φ_lo, φ_hi), (z_lo, z_hi)]\n\n\n\n\n\n","category":"function"},{"location":"api/#GravitationalPotentials.bounds_t","page":"API reference","title":"GravitationalPotentials.bounds_t","text":"bounds_t(model)\n\nCylindrical coordinates bounds ([s_lo, φ_lo, z_lo], [s_hi, φ_hi, z_hi])\n\n\n\n\n\n","category":"function"},{"location":"api/#Density-models","page":"API reference","title":"Density models","text":"","category":"section"},{"location":"api/","page":"API reference","title":"API reference","text":"MassDensityModel\nmass_density","category":"page"},{"location":"api/#GravitationalPotentials.MassDensityModel","page":"API reference","title":"GravitationalPotentials.MassDensityModel","text":"abstract type MassDensityModel\n\nAbstract type for all mass density models. Any subtype of MassDensityModel should have a corresponding mass_density implementation.\n\n\n\n\n\n","category":"type"},{"location":"api/#GravitationalPotentials.mass_density","page":"API reference","title":"GravitationalPotentials.mass_density","text":"mass_density(model::MassDensityModel, s, φ, z)\n\nGiven a mass density profile/model, return the corresponding mass density at cylindrical coordinates (s φ z).\n\n\n\n\n\n","category":"function"},{"location":"api/#Uniform-sphere-model","page":"API reference","title":"Uniform sphere model","text":"","category":"section"},{"location":"api/","page":"API reference","title":"API reference","text":"UniformSphereDensity\nmass_density(::UniformSphereDensity, ::Any, ::Any, ::Any)","category":"page"},{"location":"api/#GravitationalPotentials.UniformSphereDensity","page":"API reference","title":"GravitationalPotentials.UniformSphereDensity","text":"UniformSphereDensity <: MassDensityModel\n\nModel for a uniform density sphere centered at the origin.\n\nFields\n\nrₛ::Any: radius of the sphere\nρₛ::Any: density of the sphere\n\nSee also the implementation of mass_density for this type.\n\n\n\n\n\n","category":"type"},{"location":"api/#GravitationalPotentials.mass_density-Tuple{UniformSphereDensity, Any, Any, Any}","page":"API reference","title":"GravitationalPotentials.mass_density","text":"mass_density(model::UniformSphereDensity, s, φ, z)\n\nDensity model\n\nρ(s φ z) = begincases\n    ρ_s  textif  sqrts^2 + z^2  R_s \n    0  textotherwise\nendcases\n\nwhere\n\nρ_s = model.ρₛ\nR_s = model.rₛ\n\n\n\n\n\n","category":"method"},{"location":"api/#Power-law-sphere-model","page":"API reference","title":"Power law sphere model","text":"","category":"section"},{"location":"api/","page":"API reference","title":"API reference","text":"PowerLawSphereDensity\nmass_density(::PowerLawSphereDensity, ::Any, ::Any, ::Any)","category":"page"},{"location":"api/#GravitationalPotentials.PowerLawSphereDensity","page":"API reference","title":"GravitationalPotentials.PowerLawSphereDensity","text":"PowerLawSphereDensity <: MassDensityModel\n\nModel for a sphere with a power law density centered at the origin.\n\nFields\n\nrₛ::Any: radius of the sphere\nρ₀::Any: scale density of the sphere, density changes by r₀^α at radius = r₀\nr₀::Any: scale radius of the sphere\nα::Float64: power by which density reduces with respect to radius\n\nSee also the implementation of mass_density for this type.\n\n\n\n\n\n","category":"type"},{"location":"api/#GravitationalPotentials.mass_density-Tuple{PowerLawSphereDensity, Any, Any, Any}","page":"API reference","title":"GravitationalPotentials.mass_density","text":"mass_density(model::PowerLawSphereDensity, s, φ, z)\n\nDensity model\n\nρ(s φ z) = begincases\n    ρ_0 left(sqrts^2 + z^2R_0right)^α  textif  sqrts^2 + z^2  R_s \n    0  textotherwise\nendcases\n\nwhere\n\nρ_0 = model.ρ₀\nR_0 = model.r₀\nα = model.α\nR_s = model.rₛ\n\n\n\n\n\n","category":"method"},{"location":"api/#Uniform-cylinder-density-model","page":"API reference","title":"Uniform cylinder density model","text":"","category":"section"},{"location":"api/","page":"API reference","title":"API reference","text":"UniformCylinderDensity\nmass_density(::UniformCylinderDensity, ::Any, ::Any, ::Any)","category":"page"},{"location":"api/#GravitationalPotentials.UniformCylinderDensity","page":"API reference","title":"GravitationalPotentials.UniformCylinderDensity","text":"UniformCylinderDensity <: MassDensityModel\n\nModel for a uniform density cylinder centered at the origin.\n\nFields\n\nr_c::Any: radius of the cylinder\nh_c::Any: height of the cylinder\nρ_c::Any: density of the cylinder\n\nSee also the implementation of mass_density for this type.\n\n\n\n\n\n","category":"type"},{"location":"api/#GravitationalPotentials.mass_density-Tuple{UniformCylinderDensity, Any, Any, Any}","page":"API reference","title":"GravitationalPotentials.mass_density","text":"mass_density(model::UniformCylinderDensity, s, φ, z)\n\nDensity model\n\nρ(s φ z) = begincases\n    ρ_0  textif  s  R_c z  H_c \n    0  textotherwise\nendcases\n\nwhere\n\nρ_0 = model.ρ_c\nR_c = model.r_c\nH_c = model.h_c\n\n\n\n\n\n","category":"method"},{"location":"api/#Power-law-cylinder-density-model","page":"API reference","title":"Power law cylinder density model","text":"","category":"section"},{"location":"api/","page":"API reference","title":"API reference","text":"PowerLawCylinderDensity\nmass_density(::PowerLawCylinderDensity, ::Any, ::Any, ::Any)","category":"page"},{"location":"api/#GravitationalPotentials.PowerLawCylinderDensity","page":"API reference","title":"GravitationalPotentials.PowerLawCylinderDensity","text":"PowerLawCylinderDensity <: MassDensityModel\n\nModel for a cylinder with a power law density centered at the origin.\n\nFields\n\nr_c::Any: radius of the cylinder\nh_c::Any: height of the cylinder\nρ₀::Any: scale density of the sphere\nr₀::Any: scale radius of the cylinder\nh₀::Any: scale height of the cylinder\nα::Float64: power by which density changes with respect to radius\nβ::Float64: power by which density changes with respect to height\n\nSee also the implementation of mass_density for this type.\n\n\n\n\n\n","category":"type"},{"location":"api/#GravitationalPotentials.mass_density-Tuple{PowerLawCylinderDensity, Any, Any, Any}","page":"API reference","title":"GravitationalPotentials.mass_density","text":"mass_density(model::PowerLawCylinderDensity, s, φ, z)\n\nDensity model\n\nρ(s φ z) = begincases\n    ρ_0 left(sR_0right)^α left(zh_0right)^β  textif  s  R_c  z  H_c \n    0  textotherwise\nendcases\n\nwhere\n\nρ_0 = model.ρ₀\nR_0 = model.r₀\nH_0 = model.h₀\nα = model.α\nβ = model.β\nR_c = model.r_c\nH_c = model.h_c\n\n\n\n\n\n","category":"method"},{"location":"api/#Spiral-galaxy-model","page":"API reference","title":"Spiral galaxy model","text":"","category":"section"},{"location":"api/","page":"API reference","title":"API reference","text":"SpiralGalaxyDensity\nmass_density(::SpiralGalaxyDensity, ::Any, ::Any, ::Any)","category":"page"},{"location":"api/#GravitationalPotentials.SpiralGalaxyDensity","page":"API reference","title":"GravitationalPotentials.SpiralGalaxyDensity","text":"SpiralGalaxyDensity <: MassDensityModel\n\nModel for a spiral galaxy that has both a disk and a central bulge, with the bulge density dominating in the overlap.\n\nThe bulge is modeled as a uniform density sphere, and the disk is modeled as a uniform density cylinder.\n\nFields\n\nr_bulge::Any: radius of the bulge\nr_disk::Any: radius of the disk\nh_disk::Any: height of the disk\nρ_bulge::Any: density of the bulge\nρ_disk::Any: density of the disk\n\nSee also the implementation of mass_density for this type.\n\n\n\n\n\n","category":"type"},{"location":"api/#GravitationalPotentials.mass_density-Tuple{SpiralGalaxyDensity, Any, Any, Any}","page":"API reference","title":"GravitationalPotentials.mass_density","text":"mass_density(model::SpiralGalaxyDensity, s, φ, z)\n\nDensity model\n\nρ(s φ z) = begincases\n    ρ_b  textif  sqrts^2 + z^2  R_b \n    ρ_d  textif  s  R_d z  H_d \n    0  textotherwise\nendcases\n\nwhere\n\nρ_b = model.ρ_bulge\nρ_d = model.ρ_disk\nR_b = model.r_bulge\nR_d = model.r_disk\nH_d = model.h_disk\n\n\n\n\n\n","category":"method"},{"location":"api/#Einasto-density-model","page":"API reference","title":"Einasto density model","text":"","category":"section"},{"location":"api/","page":"API reference","title":"API reference","text":"EinastoDensity\nmass_density(::EinastoDensity, ::Any, ::Any, ::Any)","category":"page"},{"location":"api/#GravitationalPotentials.EinastoDensity","page":"API reference","title":"GravitationalPotentials.EinastoDensity","text":"EinastoDensity <: MassDensityModel\n\nFields\n\nρ₀::Any: central density\na₀::Any: harmonic mean radius\nN::Float64: structural parameter\nk::Float64: normalizing constant\n\nSee Einasto, J. and Haud, U., “Galactic models with massive corona. I. Method”, Astronomy and Astrophysics, vol. 223, no. 1, pp. 89–94, 1989.\n\nSee also the implementation of mass_density for this type.\n\n\n\n\n\n","category":"type"},{"location":"api/#GravitationalPotentials.mass_density-Tuple{EinastoDensity, Any, Any, Any}","page":"API reference","title":"GravitationalPotentials.mass_density","text":"mass_density(model::EinastoDensity, s, φ, z)\n\nDensity model\n\nρ(s φ z) = ρ₀ expleft-left(fracak a₀right)^1Nright\n\nwhere\n\na = sqrts² + z²\nρ₀ = model.ρ₀\na₀ = model.a₀\nk = model.k\nN = model.N\n\n\n\n\n\n","category":"method"},{"location":"api/#Navarro–Frenk–White-Profile","page":"API reference","title":"Navarro–Frenk–White Profile","text":"","category":"section"},{"location":"api/","page":"API reference","title":"API reference","text":"NFWDensity\nmass_density(::NFWDensity, ::Any, ::Any, ::Any)","category":"page"},{"location":"api/#GravitationalPotentials.NFWDensity","page":"API reference","title":"GravitationalPotentials.NFWDensity","text":"NFWDensity <: MassDensityModel\n\nFields\n\nρ₀::Any: Scale density\nr_s::Any: Scale radius\n\n\n\n\n\n","category":"type"},{"location":"api/#GravitationalPotentials.mass_density-Tuple{NFWDensity, Any, Any, Any}","page":"API reference","title":"GravitationalPotentials.mass_density","text":"mass_density(model::NFWDensity, s, φ, z)\n\nDensity model\n\nρ(s φ z) = fracρ₀fracrr_s left(1 + fracrr_sright)^2\n\nwhere r = sqrts^2 + z^2\n\n\n\n\n\n","category":"method"},{"location":"","page":"Home","title":"Home","text":"CurrentModule = GravitationalPotentials","category":"page"},{"location":"#GravitationalPotentials","page":"Home","title":"GravitationalPotentials","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for GravitationalPotentials.","category":"page"}]
}
