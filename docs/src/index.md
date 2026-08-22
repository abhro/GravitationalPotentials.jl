```@meta
CurrentModule = GravitationalPotentials
```

# GravitationalPotentials

Documentation for [GravitationalPotentials](https://github.com/abhro/GravitationalPotentials.jl).

This project provides mass distributions and their gravitational potential fields in Newtonian dynamics.

## Theory basics

For a mass distribution ``ρ(\mathbf{r})`` in a volume ``Ω``, the gravitational potential ``Φ`` exerted on a point
``\mathbf{r}`` by the distribution is

```math
Φ(\mathbf{r}) = \int_Ω \frac{ρ(\mathbf{r}') \, d^3\mathbf{r}'}{\left|\mathbf{r} - \mathbf{r}'\right|}
```

and the gravitational acceleration field is

```math
\mathbf{g}(\mathbf{r}) = - \boldsymbol{∇} Φ(\mathbf{r})
= \int_Ω \frac{\mathbf{r} - \mathbf{r}'}{\left|\mathbf{r} - \mathbf{r}'\right|^3} \, ρ(\mathbf{r}') \, d^3\mathbf{r}'.
```
