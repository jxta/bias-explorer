# bias_explorer — LMFDB label → Chebyshev-bias graph

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/jxta/bias-explorer/main?labpath=bias_explorer.ipynb)

> ⬆️ **Launch on Binder** to draw a bias graph from just an LMFDB label — no install needed.

Give an LMFDB elliptic-curve label (e.g. `11.a1`, `37.a1`, `389.a1`, `5077.a1`)
and get its Chebyshev-bias plot

  S(x) = Σ_{p≤x} a_p/p  ~  (1/2 − rank)·loglog x + c    (Aoki–Koyama bias).

The red dashed line is the Aoki–Koyama theory slope `1/2 − (analytic rank)`.

## Run
- **On Binder / NII Jupyter:** click the badge above. After the (one-time) build,
  open `bias_explorer.ipynb`, edit the label, run.
- **Locally:** `conda env create -f environment.yml && conda activate bias-explorer`,
  then `python bias_explorer.py 389.a1` or open `bias_explorer.ipynb`.

## Notes
- a_p computed with PARI via cypari2 (conda-forge, PARI bundled — no system install).
- Curve data fetched from the LMFDB API by label. Default x ≤ 10⁶.
- EC over Q for now; number-field / Artin versions to follow.
