#!/usr/bin/env python3
"""
bias_explorer — draw the Chebyshev-bias graph from just an LMFDB label.
(Aoki meeting / India-talk public prep: "give an LMFDB label, get the bias plot".)

Elliptic curve E/Q:  S(x) = sum_{p<=x} a_p/p  ~  (1/2 - rank)*loglog x + c
                     (Aoki-Koyama bias; motivic weight 1).

Usage:  python bias_explorer.py 11.a1
        notebook:  plot_bias("389.a1")              # default x <= 1e6
                   plot_bias("389.a1", xmax=10**8)  # larger x (slower; up to ~1e8)
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


def _sieve(n):
    """Primes up to n via a simple numpy sieve (memory ~ n bytes; fine up to ~1e8)."""
    s = np.ones(n + 1, dtype=bool); s[:2] = False
    for i in range(2, int(n**0.5) + 1):
        if s[i]:
            s[i*i::i] = False
    return np.nonzero(s)[0]


def bias_trajectory(ainvs, xmax=10**6, npts=1500):
    """S(x) = sum_{p<=x, p good} a_p/p, recorded at npts log-spaced x-thresholds
    (value at the first prime that reaches each threshold)."""
    xmax = int(xmax)
    if xmax > 2 * 10**8:
        raise ValueError("xmax too large for this in-browser demo (use <= 2e8). "
                         "Large-x runs are done on the compute nodes.")
    E = PARI.ellinit(ainvs)
    N = int(PARI.ellglobalred(E)[0])
    primes = _sieve(xmax)
    cps = np.unique(np.geomspace(20, xmax, npts).astype(np.int64))
    xs, Ss, S, ci = [], [], 0.0, 0
    for p in primes:
        p = int(p)
        if N % p:
            S += int(PARI.ellap(E, p)) / p
        while ci < len(cps) and p >= cps[ci]:
            xs.append(p); Ss.append(S); ci += 1
    return np.array(xs), np.array(Ss)


def plot_bias(label, xmax=10**6, ax=None, npts=1500):
    """LMFDB label -> Chebyshev-bias plot with the Aoki-Koyama theory line."""
    ainvs, rank, lab = fetch_ec(label)
    xs, Ss = bias_trajectory(ainvs, xmax, npts)
    if ax is None:
        _, ax = plt.subplots(figsize=(9, 5.5))
    ax.plot(xs, Ss, lw=0.6, color="C0", label=f"S(x)=$\\Sigma\\,a_p/p$  (EC {lab})")
    tail = slice(-len(Ss) // 3, None)
    c = Ss[tail].mean() - (0.5 - rank) * np.log(np.log(xs[tail])).mean()
    xg = np.geomspace(xs[0], xs[-1], 200)
    ax.plot(xg, (0.5 - rank) * np.log(np.log(xg)) + c, "r--", lw=1.2,
            label=f"theory $(1/2-\\mathrm{{rank}})\\,\\log\\log x$, rank={rank}")
    ax.set_xscale("log"); ax.set_xlabel("x"); ax.set_ylabel("S(x)")
    ax.set_title(f"Chebyshev bias of {lab}  (analytic rank {rank}),  x ≤ {xmax:g}")
    ax.grid(alpha=0.3); ax.legend()
    return ax, rank


if __name__ == "__main__":
    label = sys.argv[1] if len(sys.argv) > 1 else "11.a1"
    xmax = int(float(sys.argv[2])) if len(sys.argv) > 2 else 10**6
    ax, rank = plot_bias(label, xmax=xmax)
    out = f"/tmp/bias_{label.replace('.', '_')}.png"
    ax.figure.savefig(out, dpi=130, bbox_inches="tight")
    print(f"OK label={label} xmax={xmax} rank={rank} -> {out}")
