#!/usr/bin/env python3
"""
bias_explorer — draw the Chebyshev-bias graph from just an LMFDB label.
(Aoki meeting / India-talk public prep: "give an LMFDB label, get the bias plot".)

Elliptic curve E/Q:  S(x) = sum_{p<=x} a_p/p  ~  (1/2 - rank)*loglog x + c
                     (Aoki-Koyama bias; motivic weight 1).
Dependencies are pip-only (cypari2 = PARI, matplotlib, numpy) so it runs on
mybinder.org / NII Jupyter with no system PARI/sage needed.

Usage:  python bias_explorer.py 11.a1
        notebook:  plot_bias("389.a1")
"""
import sys, json, urllib.request
import numpy as np
import matplotlib.pyplot as plt
import cypari2

PARI = cypari2.Pari()


def fetch_ec(label):
    """Fetch (ainvs, rank, canonical_label) for an EC from the LMFDB API."""
    if "." not in label:                       # tolerant: 11a1 -> 11.a1
        i = next(k for k, c in enumerate(label) if c.isalpha())
        label = label[:i] + "." + label[i:]
    url = f"https://www.lmfdb.org/api/ec_curvedata/?lmfdb_label={label}&_format=json"
    with urllib.request.urlopen(url, timeout=30) as r:
        rec = json.load(r)["data"][0]
    return list(rec["ainvs"]), int(rec["rank"]), rec["lmfdb_label"]


def bias_trajectory(ainvs, xmax=10**6, npts=600):
    """S(x) = sum_{p<=x, p good} a_p/p, sampled at npts log-spaced checkpoints."""
    E = PARI.ellinit(ainvs)
    N = int(PARI.ellglobalred(E)[0])
    primes = [int(p) for p in PARI.primes(int(PARI.primepi(xmax)))]
    checkpoints = set(np.unique(np.geomspace(30, xmax, npts).astype(np.int64)).tolist())
    xs, Ss, S = [], [], 0.0
    for p in primes:
        if N % p:
            S += int(PARI.ellap(E, p)) / p
        if p in checkpoints:
            xs.append(p); Ss.append(S)
    return np.array(xs), np.array(Ss)


def plot_bias(label, xmax=10**6, ax=None):
    """LMFDB label -> Chebyshev-bias plot with the Aoki-Koyama theory line."""
    ainvs, rank, lab = fetch_ec(label)
    xs, Ss = bias_trajectory(ainvs, xmax)
    if ax is None:
        _, ax = plt.subplots(figsize=(9, 5.5))
    ax.plot(xs, Ss, lw=0.6, color="C0", label=f"S(x)=$\\Sigma\\,a_p/p$  (EC {lab})")
    tail = slice(-len(Ss) // 3, None)
    c = Ss[tail].mean() - (0.5 - rank) * np.log(np.log(xs[tail])).mean()
    xg = np.geomspace(xs[0], xs[-1], 200)
    ax.plot(xg, (0.5 - rank) * np.log(np.log(xg)) + c, "r--", lw=1.2,
            label=f"theory $(1/2-\\mathrm{{rank}})\\,\\log\\log x$, rank={rank}")
    ax.set_xscale("log"); ax.set_xlabel("x"); ax.set_ylabel("S(x)")
    ax.set_title(f"Chebyshev bias of {lab}  (analytic rank {rank})")
    ax.grid(alpha=0.3); ax.legend()
    return ax, rank


if __name__ == "__main__":
    label = sys.argv[1] if len(sys.argv) > 1 else "11.a1"
    ax, rank = plot_bias(label)
    out = f"/tmp/bias_{label.replace('.', '_')}.png"
    ax.figure.savefig(out, dpi=130, bbox_inches="tight")
    print(f"OK label={label} rank={rank} -> {out}")
