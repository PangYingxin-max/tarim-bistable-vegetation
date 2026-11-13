"""Structural versus observational identifiability of the two thresholds.

The two saddle-node points of

    dP/dt = rho P (P - A)(K - P),      A = A0 + A1 T - A2 R

fall where the Allee threshold crosses the two boundary branches:

    A(T*, R)  = 0   ->  R_rec = (A0 + A1 T*)       / A2 ,  T*  = s0/mu
    A(T**, R) = K   ->  R_deg = (A0 + A1 T** - K)  / A2 ,  T** solves
                        s0 - aK/(1+g T) + eta K - mu T = 0

so the bistable window is

    dR = R_rec - R_deg = (K - A1 (T** - T*)) / A2 .

Three consequences are tested here numerically:

  (1) the fold exists if and only if K > A1 (T** - T*), a condition that
      involves neither A2 nor rho, the two parameters the record cannot see;
  (2) the window width scales exactly as 1/A2, i.e. the ratio R_rec/R_deg is
      invariant to A2 while the absolute positions are not;
  (3) a Monte-Carlo perturbation of all ten parameters (+-10 %, +-20 %)
      leaves the fold intact in the great majority of draws, so the missing
      identification is observational rather than structural.
"""
import json
import sys

import numpy as np

import paths
from model_v8 import load_parameters, scan_thresholds

BASE = load_parameters()


def stress_roots(P, p, lo=1e-6, hi=3.0):
    """Roots of the stress balance G(P,T)=0 on (lo,hi), guarded for the MC."""
    def G(T):
        return (p["s0"] - p["a"] * P / (1.0 + p["g"] * T)
                + p["eta"] * P - p["mu"] * T)
    if G(lo) * G(hi) > 0:
        return []
    a, b = lo, hi
    for _ in range(200):
        m = 0.5 * (a + b)
        if G(a) * G(m) <= 0:
            b = m
        else:
            a = m
    return [0.5 * (a + b)]


def thresholds(p):
    """Closed-form saddle-node positions; empty tuple if a branch is missing."""
    tb = stress_roots(0.0, p)
    th = stress_roots(p["K"], p)
    if not tb or not th:
        return None
    Tb, Th = tb[0], th[0]
    return ((p["A0"] + p["A1"] * Th - p["K"]) / p["A2"],
            (p["A0"] + p["A1"] * Tb) / p["A2"], Tb, Th)


def fold_gap(p):
    """K - A1 (T** - T*); positive means a bistable window exists."""
    _, _, Tb, Th = thresholds(p)
    return p["K"] - p["A1"] * (Th - Tb), Tb, Th


def report(p):
    got = thresholds(p)
    if got is None:
        return None
    Rg, Rr, Tb, Th = got
    gap, _, _ = fold_gap(p)
    return dict(R_deg=Rg, R_rec=Rr, width=Rr - Rg, T_bare=Tb, T_healthy=Th,
                fold_gap=gap, ratio=(Rr / Rg if abs(Rg) > 1e-9 else float("nan")))


def mc(n=20000, frac=0.10, seed=0):
    rng = np.random.default_rng(seed)
    keys = ["rho", "K", "A0", "A1", "A2", "s0", "a", "g", "eta", "mu"]
    out = []
    for _ in range(n):
        p = dict(BASE)
        for k in keys:
            p[k] = BASE[k] * (1.0 + frac * rng.uniform(-1, 1))
        r = report(p)
        if r is not None:
            out.append(r)
    return out


if __name__ == "__main__":
    base = report(BASE)
    print("baseline")
    for k, v in base.items():
        print("  %-11s %10.4f" % (k, v))

    # (1) closed form versus numerical continuation
    found = scan_thresholds(BASE)
    print("\ncontinuation transitions (R, n_stable before, n_stable after):")
    for f in found:
        print("   %.1f  %d -> %d" % f)

    # (2) A2 scaling
    print("\nA2 scaling (A2 multiplied by f):")
    for f in (0.6, 0.8, 1.0, 1.25, 1.6):
        p = dict(BASE)
        p["A2"] = BASE["A2"] * f
        rep = report(p)
        print("   f=%.2f  R_deg=%6.2f  R_rec=%6.2f  width=%6.2f  ratio=%.3f"
              % (f, rep["R_deg"], rep["R_rec"], rep["width"], rep["ratio"]))

    # (3) fold condition: A1(T**-T*) < K
    print("\nfold condition  K - A1(T**-T*) = %.4f  (>0 means bistable)"
          % base["fold_gap"])

    # (4) Monte Carlo
    summary = {}
    for frac in (0.10, 0.20):
        reps = mc(n=12000, frac=frac, seed=int(frac * 1000))
        ok = [r for r in reps if r["width"] > 0 and r["R_deg"] > 0]
        w = np.array([r["width"] for r in ok])
        rg = np.array([r["R_deg"] for r in ok])
        rr = np.array([r["R_rec"] for r in ok])
        ra = np.array([r["ratio"] for r in ok if np.isfinite(r["ratio"])])
        summary["%.2f" % frac] = dict(
            draws=len(reps),
            bistable_fraction=len(ok) / len(reps),
            R_deg_p05=float(np.percentile(rg, 5)), R_deg_p95=float(np.percentile(rg, 95)),
            R_rec_p05=float(np.percentile(rr, 5)), R_rec_p95=float(np.percentile(rr, 95)),
            width_p05=float(np.percentile(w, 5)), width_p95=float(np.percentile(w, 95)),
            R_deg_median=float(np.median(rg)), R_rec_median=float(np.median(rr)),
            width_median=float(np.median(w)),
            ratio_p05=float(np.percentile(ra, 5)), ratio_p95=float(np.percentile(ra, 95)),
            ratio_cv=float(np.std(ra) / np.mean(ra)))
        print("\nperturbation +-%d%% over 10 parameters (%d draws)" %
              (int(frac * 100), len(reps)))
        s = summary["%.2f" % frac]
        print("   bistable window in R>0 in %.1f %% of draws"
              % (100 * s["bistable_fraction"]))
        print("   R_deg  5-95%%: %.1f - %.1f mm/yr" % (s["R_deg_p05"], s["R_deg_p95"]))
        print("   R_rec  5-95%%: %.1f - %.1f mm/yr" % (s["R_rec_p05"], s["R_rec_p95"]))
        print("   width  5-95%%: %.1f - %.1f mm/yr" % (s["width_p05"], s["width_p95"]))
        print("   medians: R_deg %.1f  R_rec %.1f  width %.1f"
              % (s["R_deg_median"], s["R_rec_median"], s["width_median"]))
        print("   ratio  5-95%%: %.2f - %.2f  (CV %.1f %%)"
              % (s["ratio_p05"], s["ratio_p95"], 100 * s["ratio_cv"]))

    out = dict(baseline=base, continuation=[list(f) for f in found],
               mc=summary,
               channel=None)
    dest = paths.OUT / "identifiability_v8.json"
    with open(dest, "w", encoding="utf-8") as fh:
        json.dump(out, fh, indent=2)
    print("\nwritten", dest)
