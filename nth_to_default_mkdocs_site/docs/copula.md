# Copulas and default dependence

We preserve marginal survival curves \(S_i(t)\) and impose dependence through a copula \(C\).

## Sklar’s theorem (construction)

Let \(F_i(t)=\mathbb{Q}(\tau_i\le t)\) be marginal default-time CDFs.
Sklar’s theorem implies

$$
\mathbb{Q}(\tau_1\le t_1,\dots,\tau_d\le t_d)
= C\!\left(F_1(t_1),\dots,F_d(t_d)\right).
$$

Simulate \(U=(U_1,\dots,U_d)\) from \(C\) and map to default times using the inverse cumulative hazard.

## Gaussian copula

Let \(Z\sim\mathcal{N}(0,\Sigma)\) with correlation matrix \(\Sigma\). Define

$$
U_i = \Phi(Z_i),
$$

where \(\Phi\) is the standard normal CDF.

## Student-t copula (tail dependence)

Scale-mixture construction:

- \(Z\sim\mathcal{N}(0,\Sigma)\)
- \(V\sim \chi^2_\nu\), \(W=\sqrt{V/\nu}\)
- \(T = Z/W\)

Then

$$
U_i = t_\nu(T_i).
$$

Smaller \(\nu\) yields stronger symmetric tail dependence.
