# Correlated Nth-to-Default Basket CDS

This site documents a market-consistent valuation framework for **kth-to-default (Nth-to-default) basket CDS** under
piecewise-constant intensities and **elliptical copulas** (Gaussian / Student-t), with an explicit bridge to **risk-capital** analytics.

!!! note "How equations render"
    This documentation is configured with **MathJax**. Inline math uses `$...$` and display math uses `$$...$$`.

## Quick start (local docs)

```bash
pip install -r requirements-docs.txt
mkdocs serve
```

## Project map

1. **Intensity calibration**: bootstrap a piecewise-constant hazard rate term structure from CDS par spreads.
2. **Copula sampling**: generate correlated uniforms $(U_1,\dots,U_d)$ via Gaussian or Student-$t$ copula.
3. **Default time inversion**: map uniforms to default times using inverse cumulative hazard.
4. **Kth trigger and cashflows**: compute $\tau_{(k)}$ and value protection/premium legs.
5. **Risk capital**: convert scenario engine into horizon loss distribution and compute EL/VaR/ES.
