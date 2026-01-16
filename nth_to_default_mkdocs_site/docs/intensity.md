# Intensity calibration from single-name CDS

This section formalizes the reduced-form marginal model and the bootstrap of a **piecewise-constant hazard rate** from CDS par spreads.

## Reduced-form setup

For obligor \(i\), let \(\tau_i\) denote default time and \(\lambda_i(t)\) the (risk-neutral) default intensity.
Define the cumulative hazard

$$
H_i(t) = \int_0^t \lambda_i(s)\,ds,
$$

so the survival function is

$$
S_i(t) = \mathbb{Q}(\tau_i > t) = \exp\{-H_i(t)\}.
$$

A standard simulation identity is: if \(U_i \sim \mathrm{Unif}(0,1)\), then

$$
\tau_i = H_i^{-1}\!\left(-\log U_i\right)
$$

has survival \(S_i(t)\).

## CDS discretization (premium and protection legs)

Let coupon dates be \(0=t_0 < t_1 < \dots < t_m=T\) with accrual fractions \(\Delta_j\).
Assume deterministic discounting \(DF(t)\) and constant recovery \(R\) (LGD \(=1-R\)).

A common discrete approximation to the protection leg is

$$
PV_{\text{prot}} \approx (1-R)\sum_{j=1}^m DF(t_j)\bigl(S(t_{j-1})-S(t_j)\bigr).
$$

The premium leg per unit spread (PV01) is

$$
PV01 \approx \sum_{j=1}^m DF(t_j)\Delta_j S(t_j) \;+\; PV_{\text{accr}},
$$

where \(PV_{\text{accr}}\) is the accrued premium on default.

The par spread \(s^\*\) at maturity \(T\) satisfies

$$
s^\* = \frac{PV_{\text{prot}}}{PV01}.
$$

## Piecewise-constant hazard bootstrap

Let market maturities be \(T_1<\dots<T_n\). Assume

$$
\lambda(t) = \lambda_k,\quad t\in(T_{k-1},T_k],\quad T_0=0.
$$

Then for \(t\in(T_{k-1},T_k]\),

$$
H(t) = \sum_{j=1}^{k-1}\lambda_j (T_j-T_{j-1}) + \lambda_k (t-T_{k-1}).
$$

Bootstrap proceeds sequentially by solving for each segment hazard \(\lambda_k\) to match the par spread at \(T_k\).

## Inversion of a piecewise hazard

Given \(x=-\log U\), locate the first segment where \(H(T_k)\ge x\). Then solve

$$
t = T_{k-1} + \frac{x - H(T_{k-1})}{\lambda_k}.
$$
