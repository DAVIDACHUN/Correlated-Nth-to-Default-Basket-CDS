# Nth-to-Default Basket CDS (Copula / Intensity) – Python Refactor

This repo contains a clean Python implementation of an Nth-to-Default basket CDS model:
- market-implied intensity curves bootstrapped from CDS spreads,
- Gaussian and Student-t copulas for default dependence,
- Monte Carlo pricing of kth-to-default fair spreads,
- (optional) simple loss distribution metrics useful for risk-capital narratives (EL / VaR / ES).

## Quick start (local)

```bash
python -m venv .venv
source .venv/bin/activate     # mac/linux
# .venv\Scripts\activate    # windows
pip install numpy pandas scipy
python -c "from basket_cds import *; print('ok')"
```

## Example usage

See `example_run.py` (included) for a full worked example using the June-2018 single-name CDS spreads from the paper.

## GitHub Pages (site hosting)

If you want a lightweight project site (docs + figures + links):

1. Create a repo on GitHub (e.g., `nth-to-default-basket-cds`).
2. Push this folder to the repo (see commands below).
3. Add a `docs/` folder with an `index.html` (or use MkDocs).
4. In GitHub:
   - Settings → Pages
   - Source: `Deploy from a branch`
   - Branch: `main`
   - Folder: `/docs`
5. Your site will publish at: `https://<username>.github.io/<repo>/`

### Minimal Pages site

Create `docs/index.html`:

```html
<!doctype html>
<html>
  <head><meta charset="utf-8"><title>Nth-to-Default Basket CDS</title></head>
  <body>
    <h1>Nth-to-Default Basket CDS</h1>
    <p>Python implementation + paper.</p>
    <ul>
      <li><a href="../nth_to_default_polished.docx">Paper (DOCX)</a></li>
      <li><a href="../basket_cds.py">Python module</a></li>
    </ul>
  </body>
</html>
```

Commit and push, then enable Pages.

## Git commands (copy/paste)

```bash
git init
git add .
git commit -m "Initial commit: Nth-to-default basket CDS (Python refactor)"
git branch -M main
git remote add origin https://github.com/<username>/<repo>.git
git push -u origin main
```

> Tip: If you also want a Jupyter notebook rendered nicely on Pages, convert it to HTML and place it in `docs/`.
