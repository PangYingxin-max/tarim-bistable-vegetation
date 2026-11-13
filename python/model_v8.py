"""Vegetation--water model calibrated to the 2000-2023 MODIS record.

    dP/dt = rho P (P - A) (K - P)
    dT/dt = s0 - a P/(1 + g T) + eta P - mu T
    A(T,R) = A0 + A1 T - A2 R

The effective water input R acts on the persistence threshold A (root-zone
water available for establishment) and not on the surface-stress balance: the
observed TVDI of the corridor has no trend while R doubled, so a model in
which R only lowers T cannot reproduce the transition.

The baseline values are the published ones (Table 1 of the paper) and are also
stored in data/calibration_v10.json.  Run this module directly to print the
two thresholds and the equilibria.
"""
import json
from pathlib import Path

import numpy as np

DATA = Path(__file__).resolve().parents[1] / "data"

# Published baseline (Table 1, Table S2).  These reproduce R_deg = 36.0 and
# R_rec = 55.0 mm/yr, the values used throughout the manuscript.
BASE = dict(rho=0.80, K=0.41, A0=-1.6370829, A1=3.00, A2=0.01893878,
            s0=0.71432416, a=1.00, g=0.50, eta=0.72, mu=0.80)


def load_parameters(path=None):
    """Baseline parameters from data/calibration_v10.json, else BASE."""
    path = Path(path) if path else DATA / "calibration_v10.json"
    if path.exists():
        with open(path, encoding="utf-8") as fh:
            return {k: float(v) for k, v in json.load(fh)["parameters"].items()}
    return dict(BASE)


def A_of(T, R, p):
    return p["A0"] + p["A1"] * T - p["A2"] * R


def F(P, T, R, p):
    return p["rho"] * P * (P - A_of(T, R, p)) * (p["K"] - P)


def G(P, T, p):
    return p["s0"] - p["a"] * P / (1 + p["g"] * T) + p["eta"] * P - p["mu"] * T


def jac(P, T, R, p):
    A, K = A_of(T, R, p), p["K"]
    j11 = p["rho"] * ((P - A) * (K - P) + P * (K - P) - P * (P - A))
    j12 = -p["rho"] * p["A1"] * P * (K - P)
    j21 = -p["a"] / (1 + p["g"] * T) + p["eta"]
    j22 = p["a"] * p["g"] * P / (1 + p["g"] * T) ** 2 - p["mu"]
    return np.array([[j11, j12], [j21, j22]])


def T_of_P(P, p, hi=1.5):
    """Stress roots of G(P,T)=0 on (0,hi) (G is monotone decreasing in T).

    T is a normalized index: the observed corridor sits at 0.90, close to its
    upper bound, so perturbed parameter sets are allowed to run past 1.
    """
    lo = 1e-6
    if G(P, lo, p) * G(P, hi, p) > 0:
        return []
    a, b = lo, hi
    for _ in range(120):
        m = 0.5 * (a + b)
        if G(P, a, p) * G(P, m, p) <= 0:
            b = m
        else:
            a = m
    return [0.5 * (a + b)]


def equilibria(R, p):
    out = []
    for P in (0.0, p["K"]):
        for T in T_of_P(P, p):
            out.append((P, T))
    # interior branch: P = A(T,R) with G(A(T,R),T)=0
    Ts = np.linspace(1e-4, 0.999, 4000)
    vals = [G(A_of(t, R, p), t, p) for t in Ts if 0 < A_of(t, R, p) < p["K"]]
    Tk = [t for t in Ts if 0 < A_of(t, R, p) < p["K"]]
    for i in range(len(vals) - 1):
        if vals[i] * vals[i + 1] < 0:
            a, b, fa = Tk[i], Tk[i + 1], vals[i]
            for _ in range(80):
                m = 0.5 * (a + b)
                fm = G(A_of(m, R, p), m, p)
                if fa * fm <= 0:
                    b = m
                else:
                    a, fa = m, fm
            out.append((A_of(0.5 * (a + b), R, p), 0.5 * (a + b)))
    uniq = []
    for P, T in out:
        if not any(abs(P - q) < 1e-6 and abs(T - s) < 1e-6 for q, s in uniq):
            uniq.append((P, T))
    res = []
    for P, T in uniq:
        ev = np.linalg.eigvals(jac(P, T, R, p))
        res.append(dict(P=P, T=T, stable=bool(np.all(ev.real < 0)), eig=ev))
    return sorted(res, key=lambda d: d["P"])


def thresholds(p):
    Tb = T_of_P(0.0, p)[0]
    Th = T_of_P(p["K"], p)[0]
    R_rec = (p["A0"] + p["A1"] * Tb) / p["A2"]          # A = 0 on the bare branch
    R_deg = (p["A0"] + p["A1"] * Th - p["K"]) / p["A2"]  # A = K on the healthy branch
    return R_deg, R_rec, Tb, Th


def scan_thresholds(p, lo=0.0, hi=100.0, n=2001):
    """Numerical (continuation) version, independent of the closed form."""
    prev = None
    found = []
    for R in np.linspace(lo, hi, n):
        eq = equilibria(R, p)
        n_stable = sum(e["stable"] for e in eq)
        if prev is not None and n_stable != prev:
            found.append((R, prev, n_stable))
        prev = n_stable
    return found


if __name__ == "__main__":
    p = BASE
    Rg, Rr, Tb, Th = thresholds(p)
    print("T_bare   = %.4f" % Tb)
    print("T_healthy= %.4f" % Th)
    print("R_deg    = %.2f mm/yr" % Rg)
    print("R_rec    = %.2f mm/yr" % Rr)
    print("window   = %.2f mm/yr" % (Rr - Rg))
    for R in (28.0, 33.1, 41.0, 52.0, 65.1):
        print("\nR = %.1f" % R)
        for e in equilibria(R, p):
            print("   P=%.4f  T=%.4f  %s" % (e["P"], e["T"], "stable" if e["stable"] else "UNSTABLE"))
