# bias_explorer — LMFDB label → Chebyshev-bias graph

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/jxta/bias-explorer/main?labpath=bias_explorer.ipynb)

> ⬆️ **Launch on Binder** to draw a bias graph from just an LMFDB label — no install needed.

Give an LMFDB elliptic-curve label (e.g. `11.a1`, `37.a1`, `389.a1`, `5077.a1`)
and get its Chebyshev-bias plot

  S(x) = Σ_{p≤x} a_p/p  ~  (1/2 − rank)·loglog x + c    (Aoki–Koyama bias).

The red dashed line is the Aoki–Koyama theory slope `1/2 − (analytic rank)`.

## Run
- **On Binder / NII Jupyter:** click the badge above (or open `bias_explorer.ipynb`),
  edit the label, run. Dependencies are pip-only via `requirements.txt` + `apt.txt`.
- **Locally:** `pip install -r requirements.txt` (needs PARI: `apt install pari-gp libpari-dev`),
  then `python bias_explorer.py 389.a1` or open `bias_explorer.ipynb`.

## Notes
- a_p computed with PARI (`ellap`); curve data fetched from the LMFDB API by label.
- Default x ≤ 10⁶. Increase `xmax` for sharper convergence.
- EC over Q for now; number-field / Artin versions to follow.
