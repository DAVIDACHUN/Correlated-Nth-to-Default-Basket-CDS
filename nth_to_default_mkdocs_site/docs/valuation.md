# Kth-to-default valuation

Let \(d\) be the number of names and \(\tau_{(k)}\) the \(k\)-th order statistic of default times.

## Protection leg

$$
\Pi_{\text{prot}} = (1-R)\frac{N}{d}\, DF(\tau_{(k)}) \mathbf{1}\{\tau_{(k)}\le T\}.
$$

## Premium leg (PV01)

$$
PV01(\omega) = \sum_{j:t_j\le \min(\tau_{(k)},T)} DF(t_j)\Delta_j \;+\; PV_{\text{accr}}(\omega).
$$

Pathwise accrued premium if \(\tau_{(k)}\le T\):

$$
PV_{\text{accr}}(\omega) = DF(\tau_{(k)})(\tau_{(k)}-t_{j^\*}).
$$

## Fair spread

$$
s^\* = \frac{PV_{\text{prot}}}{PV01}.
$$
