# Azimuthal integration

For an azimuthally symmetric mass distribution $ρ(s, z)$, the potential will
have a factor

```math
I = ∫_0^{2π} \frac{dφ'}{|\mathbf{r} - \mathbf{r}'|}
= ∫_0^{2π} \frac{dφ'}{\sqrt{s^2 + s'^2 - 2 s s' \cos(φ - φ') + (z - z')^2}}
```
where ``\mathbf{r} = (s, φ, z)`` is the test (observation) point and
``\mathbf{r}' = (s', φ', z')`` are the source points being integrated over.

Since the mass distribution is azimuthally symmetric, we can rotate our
coordinates to always have the test point at ``φ = 0``, so we can simplify the
integral as

```math
I = ∫_0^{2π} \frac{dφ'}{\sqrt{k^2 - p \cos(φ')}}
```

where ``k^2 = s^2 + s'^2 + (z - z')^2`` and ``p = 2 s s'``.

Here, we will attempt to express ``I`` as an incomplete elliptic integral of the first
kind (the following expression is from Wolfram|Alpha):

```math
\begin{align*}
I &= \frac{2}{\sqrt{k^2 - p}} \left.F{\left(\frac{φ'}{2}, -\frac{2p}{k^2-p}\right)}\right|_0^{2π} \\
  &= \frac{2}{\sqrt{k^2 - p}} \left[
      F{\left(π, -\frac{2p}{k^2-p}\right)} - F{\left(0, -\frac{2p}{k^2-p}\right)}
  \right]
\end{align*}
```

or maybe a complete elliptic integral too. The definitions of the integrals are

```math
\begin{align*}
F(φ, a) &= ∫_0^φ \frac{dφ'}{\sqrt{1 - a^2 \sin^2(φ')}} \\
K(a) &= F{\left(\tfrac{π}{2}, a\right)} - F{\left(0, a\right)}
= F{\left(\tfrac{π}{2}, a\right)}
= ∫_0^{π/2} \frac{dφ'}{\sqrt{1 - a^2 \sin^2(φ')}} \\
\end{align*}
```
