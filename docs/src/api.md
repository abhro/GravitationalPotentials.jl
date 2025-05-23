
# API reference

Documentation for [GravitationalPotentials](https://github.com/abhro/GravitationalPotentials.jl).

```@docs
bounds
bounds_t
```

## Density models

```@docs
MassDensityModel
mass_density
```

### Uniform sphere model
```@docs
UniformSphereDensity
mass_density(::UniformSphereDensity, ::Any, ::Any, ::Any)
```

### Power law sphere model
```@docs
PowerLawSphereDensity
mass_density(::PowerLawSphereDensity, ::Any, ::Any, ::Any)
```

### Uniform cylinder density model
```@docs
UniformCylinderDensity
mass_density(::UniformCylinderDensity, ::Any, ::Any, ::Any)
```

### Power law cylinder density model
```@docs
PowerLawCylinderDensity
mass_density(::PowerLawCylinderDensity, ::Any, ::Any, ::Any)
```

### Spiral galaxy model
```@docs
SpiralGalaxyDensity
mass_density(::SpiralGalaxyDensity, ::Any, ::Any, ::Any)
```

### Einasto density model

```@docs
EinastoDensity
mass_density(::EinastoDensity, ::Any, ::Any, ::Any)
```

### Navarro–Frenk–White Profile
```@docs
NFWDensity
mass_density(::NFWDensity, ::Any, ::Any, ::Any)
```
