# Potential due to an uniform density cylinder

## Case 1: Test point on-axis

We use the cylindrical coordinate system, where ``\mathbf{r} = (s, φ, z)`` (as in Griffiths' Introduction to
Electrodynamics). By on-axis, we mean that our test point has the constraint ``s = 0``, or, ``\mathbf{r} = (0, φ, z)``.
Note that technically ``φ`` is not well-defined for ``s = 0``, but since we have azimuthal symmetry in this system, this
hiccup will be irrelevant.

In general, we can get the gravitational potential as
```math
\Phi(\mathbf{r}) = -G \int_\text{cylinder} \frac{ρ(\mathbf{r}')}{\left|\mathbf{r} - \mathbf{r}'\right|} \, d^3 \mathbf{r}'
```
where each ``\mathbf{r}' = (s', φ', z')`` is a source point inside the cylinder, and the density function is
```math
ρ(\mathbf{r}) = \begin{cases}
ρ_0, & \text{if } \mathbf{r} ∈ \text{cylinder} \\
0,   & \text{if } \mathbf{r} ∉ \text{cylinder}
\end{cases}
```

Also note, in the cylindrical
coordinate system, we have the two relations:
```math
\begin{align*}
\left|\mathbf{r} - \mathbf{r}'\right| &= \sqrt{s'^2 + (z - z')^2}, &
d^3 \mathbf{r}' &= s' \, ds' \, dφ' \, dz'
\end{align*}
```

If we characterize the cylinder as having radius ``R`` and height ``2H``, the integral becomes
```math
\begin{align*}
\Phi(\mathbf{r}) &= - G \int_0^{2π} dφ' \int_0^R s' \, ds' \int_{-H}^{H} dz' \frac{ρ_0}{\sqrt{s'^2 + (z - z')^2}} \\[1ex]
&= - 2π G ρ_0 \int_{-H}^H dz' \int_0^R \frac{s' \, ds'}{\sqrt{s'^2 + (z - z')^2}} \\[1ex]
&= - 2π G ρ_0 \int_{-H}^H dz' {\left[\sqrt{s'^2 + (z - z')^2}\right]}_{s'=0}^{s'=R}
\end{align*}
```
where we made use of the tabulated integral
```math
\int \frac{x \, dx}{\sqrt{x^2 + k}} = \sqrt{x^2 + k} + C.
```

Plugging in the bounds for ``s'``, we have
```math
\Phi(\mathbf{r}) = -2π G ρ_0 \int_{-H}^H dz' \left(\sqrt{R^2 + (z - z')^2} - \left|z - z'\right|\right)
```

We can now use two more tabulated integrals:
```math
\begin{align*}
\int \sqrt{k + (p - x)^2} dx
    &= \frac{x-p}{2} \sqrt{k + (p - x)^2} - \frac{k^2}{2} \operatorname{artanh} \left(\frac{p-x}{\sqrt{k + (p - x)^2}}\right) + C \\
\int \left|p-x\right| dx &= \frac{x (2p-x) \operatorname{sgn}(p-x)}{2} + C
\end{align*}
```
and we get
```math
\begin{align*}
\Phi(\mathbf{r}) = - 2π G ρ_0 \Biggl[
    &\frac{z' - z}{2} \sqrt{R^2 + (z-z')^2} \\
    &- \frac{R^2}{2} \operatorname{artanh} \left(\frac{z-z'}{\sqrt{R^2 + (z-z')^2}}\right)
    - \frac{z' (2z - z') \operatorname{sgn}(p-x)}{2}
\Biggr]_{z'=-H}^{z'=H}
\end{align*}
```
```math
\begin{align*}
\Phi(\mathbf{r}) = -G π ρ_0 &\Biggl[
\hphantom{-}
\hphantom{-}
(H - z) \sqrt{R^2 + (z-H)^2} - R^2 \operatorname{artanh} \frac{z-H}{\sqrt{R^2 + (z-H)^2}} - H (2z-H) \operatorname{sgn}(z-H)
\\
&
- (-H-z) \sqrt{R^2 + (z+H)^2} + R^2 \operatorname{artanh} \frac{z+H}{\sqrt{R^2 + (z+H)^2}} - H (2z+H) \operatorname{sgn}(z+H)
\Biggr]
\end{align*}
```
```math
\begin{align*}
\Phi(\mathbf{r}) = - G π ρ_0 & \Biggl[
    (z+H) \sqrt{R^2 + (z+H)^2} - (z-H) \sqrt{R^2 + (z-H)^2}
\\ &
    + R^2 \left(\operatorname{artanh}\frac{z+H}{\sqrt{R^2 + (z+H)^2}} - \operatorname{artanh}\frac{z-H}{\sqrt{R^2+(z-H)^2}}\right)
\\ &
    - H \left((2z+H) \operatorname{sgn}(z+H) + (2z-H) \operatorname{sgn}(z-H)\right)
\Biggr]
\end{align*}
```

## Case 2: test-point on the xy-plane

Consider the point ``\mathbf{r} = (s, φ, 0)``.

```math
\begin{align*}
\Phi(\mathbf{r}) = -G ∫_{-H}^H dz' ∫_0^R ds' ∫_0^{2π} s' \, dφ' \frac{ρ_0}{\sqrt{s^2 + s'^2 - 2 s s' \cos(φ - φ') + z'^2}}
\end{align*}
```

Tabulated integral
```math
∫ \frac{dφ'}{\sqrt{K^2 - 2 s s' \cos(φ - φ')}} = - \frac{2}{\sqrt{K^2 - 2 s s'}} F{\left(\frac{φ-φ'}{2} \bigg| \frac{4 s s'}{2 s s' - K^2}\right)}
```
where ``F`` is the hypergeometric function.

And
```math
\begin{align*}
∫_{-H}^H \frac{dz'}{\sqrt{A^2 + z'^2}}
&= \left.\operatorname{artanh}\left(\frac{z'}{\sqrt{A^2 + z'^2}}\right)\right|_{-H}^H \\
&= \operatorname{artanh}\left(\frac{H}{\sqrt{A^2 + H^2}}\right)
 - \operatorname{artanh}\left(-\frac{H}{\sqrt{A^2 + H^2}}\right) \\
&= 2 \operatorname{artanh}\left(\frac{H}{\sqrt{A^2 + H^2}}\right)
\end{align*}
```
where the last equality comes from the fact that ``\operatorname{artanh}`` is an odd function.

TODO: rest of the derivation.
