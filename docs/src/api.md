```@meta
CurrentModule = GravitationalPotentials
```
# API reference

Documentation for [GravitationalPotentials](https://github.com/abhro/GravitationalPotentials.jl).

## Density models

```@docs
MassDensityModel
mass_density
```

### Uniform sphere model
```@docs
UniformSphereDensity
mass_density(::UniformSphereDensity{L}, ::L, ::Any, ::L) where {L}
```

### Power law sphere model
```@docs
PowerLawSphereDensity
mass_density(::PowerLawSphereDensity{L}, ::L, ::Any, ::L) where {L}
```

### Uniform cylinder density model
```@docs
UniformCylinderDensity
mass_density(::UniformCylinderDensity{L}, ::L, ::Any, ::L) where {L}
onaxispotential(::UniformCylinderDensity{L}, ::L) where {L}
```

### Power law cylinder density model
```@docs
PowerLawCylinderDensity
mass_density(::PowerLawCylinderDensity{L}, ::L, ::Any, ::L) where {L}
```

### Spiral galaxy model
```@docs
SpiralGalaxyDensity
mass_density(::SpiralGalaxyDensity{L}, ::L, ::Any, ::L) where {L}
```

### Einasto density model

```@docs
EinastoDensity
mass_density(::EinastoDensity{L}, ::L, ::Any, ::L) where {L}
```

### Navarro–Frenk–White Profile
```@docs
NFWDensity
mass_density(::NFWDensity{L}, ::L, ::Any, ::L) where {L}
```
