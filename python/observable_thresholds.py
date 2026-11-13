"""Translate the two thresholds into quantities a manager can observe.

Writes outputs/table_observable_thresholds.csv and prints the LaTeX rows.
"""
import io
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import model_v8 as m
import paths

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

NDVI_MIN, NDVI_MAX = 0.053, 0.230
P = m.load_parameters()
R0, KAPPA = 28.0, 5.0


def ndvi_from_P(p):
    return NDVI_MIN + (NDVI_MAX - NDVI_MIN) * p


def saddle(R):
    eq = [e for e in m.equilibria(R, P) if not e["stable"]]
    return eq[0] if eq else None


rows = []
for R in (36, 38, 40, 42, 45, 48, 50, 52, 55):
    s = saddle(R)
    pstar = s["P"] if s else 0.0
    tstar = s["T"] if s else float("nan")
    rows.append(dict(R=R, W=(R - R0) / KAPPA, P=pstar, NDVI=ndvi_from_P(pstar), T=tstar))

with open(paths.OUT / "table_observable_thresholds.csv", "w", encoding="utf-8") as fh:
    fh.write("R_mm_per_yr,W_1e8_m3_per_yr,P_unstable,NDVI_critical,T_unstable\n")
    for r in rows:
        fh.write("%.1f,%.2f,%.4f,%.4f,%.4f\n" % (r["R"], r["W"], r["P"], r["NDVI"], r["T"]))

for r in rows:
    print("%5.1f & %4.2f & %.3f & %.3f \\\\" % (r["R"], r["W"], r["P"], r["NDVI"]))

print()
print("thresholds: W_deg = %.2f, W_rec = %.2f (1e8 m3/yr)" % ((36 - R0) / KAPPA, (55 - R0) / KAPPA))
print("present R = 65.1 -> W = %.2f" % ((65.1 - R0) / KAPPA))
print("observed NDVI 2000 = 0.0758, 2023 = 0.1141")
