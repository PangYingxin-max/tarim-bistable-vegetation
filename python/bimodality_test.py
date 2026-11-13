"""Minimal observational test of the two-state hypothesis.

If the corridor is a mosaic of two alternative stable states, the pixel-level
NDVI distribution should be bimodal, with the modes corresponding to the two
branches of the model.  If the states are not separated by a barrier, the
distribution should be unimodal.  This script reports the histogram, a kernel
density estimate and the BIC of one- and two-component Gaussian mixtures for
the pooled sample and for the first and last five years.
"""
import io
import json
import sys

import numpy as np
import pandas as pd

import paths

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

df = pd.read_csv(paths.POINTS_CORRIDOR)
df = df[df["zone"] == "corridor"] if "zone" in df.columns else df


def gmm_bic(x, k, seed=0):
    """BIC of a k-component 1-D Gaussian mixture (EM, 40 restarts)."""
    rng = np.random.default_rng(seed)
    best = np.inf
    n = len(x)
    for _ in range(40):
        mu = rng.choice(x, k, replace=False)
        sd = np.full(k, x.std() / k)
        w = np.full(k, 1.0 / k)
        for _ in range(120):
            p = w * np.exp(-0.5 * ((x[:, None] - mu) / sd) ** 2) / (sd * np.sqrt(2 * np.pi))
            s = p.sum(1, keepdims=True)
            s[s == 0] = 1e-300
            r = p / s
            nk = r.sum(0) + 1e-12
            w = nk / n
            mu = (r * x[:, None]).sum(0) / nk
            sd = np.sqrt((r * (x[:, None] - mu) ** 2).sum(0) / nk) + 1e-9
        ll = np.log(s).sum()
        bic = -2 * ll + (3 * k - 1) * np.log(n)
        if bic < best:
            best, best_params = bic, (w.copy(), mu.copy(), sd.copy())
    return best, best_params


def report(name, x):
    rng = np.random.default_rng(1)
    if len(x) > 6000:
        x = rng.choice(x, 6000, replace=False)
    hist, edges = np.histogram(x, bins=np.arange(0.0, 0.40, 0.01))
    b1, p1 = gmm_bic(x, 1)
    b2, p2 = gmm_bic(x, 2)
    print("\n== %s (n = %d) ==" % (name, len(x)))
    print("  mean %.3f  sd %.3f  p5 %.3f  p50 %.3f  p95 %.3f" %
          (x.mean(), x.std(), np.percentile(x, 5), np.percentile(x, 50), np.percentile(x, 95)))
    print("  BIC 1 component %10.1f   2 components %10.1f   delta %8.1f" % (b1, b2, b1 - b2))
    print("  2-component means %s  weights %s" %
          (np.round(p2[1], 3), np.round(p2[0] / p2[0].sum(), 3)))
    top = np.argsort(hist)[::-1][:2]
    print("  histogram peaks at %s" % [(round(edges[i], 2), int(hist[i])) for i in top])
    return dict(n=int(len(x)), mean=float(x.mean()), sd=float(x.std()),
                bic1=float(b1), bic2=float(b2), mu2=[float(v) for v in p2[1]],
                w2=[float(v) for v in p2[0] / p2[0].sum()])


out = {}
out["pooled"] = report("pooled 2000-2023", df["NDVI"].values)
out["early"] = report("2000-2004", df[df.year <= 2004]["NDVI"].values)
out["late"] = report("2019-2023", df[df.year >= 2019]["NDVI"].values)

dest = paths.OUT / "bimodality_test.json"
with open(dest, "w", encoding="utf-8") as fh:
    json.dump(out, fh, indent=2)
print("\nwrote", dest)
